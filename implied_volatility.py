#!/usr/bin/env python

"""
Example for Implied Volatility using nag4py
Finds a zero of the Black Scholes function using c05auc and s30aac
Data needs to be downloaded from:
http://www.cboe.com/delayedquote/QuoteTableDownload.aspx
"""

import os, sys
import pandas
import matplotlib.pylab as plt
from nag4py.s import s30aac
from nag4py.c05 import c05auc, NAG_C05AUC_FUN
from nag4py.util import NagError, Nag_Comm, Nag_RowMajor, Nag_Call, Nag_Put, Pointer

__author__ = "John Morrissey, Chris Seymour, and Brian Spector"
__copyright__ = "Copyright 2013, The Numerical Algorithms Group Inc"
__email__ = "support@nag.co.uk"


def callback(x, comm):
    """
    Callback that calculates the Black Scholes Option Price for a given Volatility
    """

    fail = NagError()

    p_userdata = pandas.np.ctypeslib.ctypes.cast(comm[0].p, pandas.np.ctypeslib.ctypes.py_object)
    userdata = p_userdata.value

    time = pandas.np.ctypeslib.ctypes.c_double(userdata[0])
    callput = userdata[1]
    strike = pandas.np.ctypeslib.ctypes.c_double(userdata[2])
    underlying = userdata[3]
    current_price = userdata[4]
    out = pandas.np.ctypeslib.ctypes.c_double(0.0)
    
    if (x <= 0):
         x=.00001
    s30aac(Nag_RowMajor, callput, 1, 1, strike, underlying, time, x, 0.0, 0.0, out, fail)
    if(fail.code==0):
        return out.value - current_price
    #print fail.message
    return 0.0

def calcvol(exp, strike, todays_date, underlying, current_price, callput):
    """
    Root-finding method that calls NAG Library to calculate Implied Volatility
    """

    fail = NagError()

    volatility = pandas.np.ctypeslib.ctypes.c_double(.5)
    time = (exp - todays_date) / 365.0
    lb = pandas.np.ctypeslib.ctypes.c_double(0.0)
    ub = pandas.np.ctypeslib.ctypes.c_double(0.0)

    userdate = time, callput, strike, underlying, current_price
    comm = Nag_Comm()
    comm.p = pandas.np.ctypeslib.ctypes.cast(id(userdate), Pointer)

    pyfun = NAG_C05AUC_FUN(callback)

    c05auc(volatility, .1, .00001, 0.0, pyfun, lb, ub, comm, fail)

    if (fail.code == 0):
        return volatility.value
    else:
        return 0.0

# Set to hold expiration dates
dates = []

cumulative_month = {'Jan': 31, 'Feb': 57, 'Mar': 90,
                    'Apr': 120, 'May': 151, 'Jun': 181,
                    'Jul': 212, 'Aug': 243, 'Sep': 273,
                    'Oct': 304, 'Nov': 334, 'Dec': 365}

def getexpiration(x):
    monthday = x.split()
    adate = monthday[0] + ' ' + monthday[1]
    if adate not in dates:
        dates.append(adate)
    return (int(monthday[0]) - 13) * 365 + cumulative_month[monthday[1]]


def getstrike(x):
    monthday = x.split()
    return float(monthday[2])


def main():

    try:
        QuoteData = 'QuoteData.dat'
#    except IndexError:
#        sys.stderr.write("Usage: imp_vol.py QuotaData.dat\n")
#        sys.exit(1)

#        if os.path.isfile(QuoteData):
        qd = open(QuoteData, 'r')
        qd_head = []
        qd_head.append(qd.readline())
        qd_head.append(qd.readline())
        qd.close()
    except:
        sys.stderr.write("Couldn't read %s" % QuoteData)

    print "Implied Volatility for %s %s" % (qd_head[0].strip(), qd_head[1])

    # Parse the header information in QuotaData
    first = qd_head[0].split(',')
    second = qd_head[1].split()
    qd_date = qd_head[1].split(',')[0]
    
    company = first[0]
    underlyingprice = float(first[1])
    month, day = second[:2]
    today = cumulative_month[month] + int(day) - 30

    data = pandas.io.parsers.read_csv(QuoteData, sep=',', header=2, na_values=' ')

    # Need to fill the NA values in dataframe
    data = data.fillna(0.0)

    # Let's look at data where there was a recent sale 
    data = data[(data['Last Sale'] > 0) | (data['Last Sale.1'] > 0)]

    # Get the Options Expiration Date
    exp = data.Calls.apply(getexpiration)
    exp.name = 'Expiration'

    # Get the Strike Prices
    strike = data.Calls.apply(getstrike)
    strike.name = 'Strike'

    data = data.join(exp)
    data = data.join(strike)

    print 'Calculating Implied Vol of Calls...'
    impvolcall = pandas.Series(pandas.np.zeros(len(data.index)),
                               index=data.index, name='impvolCall')
    for i in data.index:
        impvolcall[i] = (calcvol(data.Expiration[i],
                                 data.Strike[i],
                                 today,
                                 underlyingprice,
                                 (data.Bid[i] + data.Ask[i]) / 2, Nag_Call))

    print 'Calculated Implied Vol for %d Calls' % len(data.index)
    data = data.join(impvolcall)

    print 'Calculating Implied Vol of Puts...'
    impvolput = pandas.Series(pandas.np.zeros(len(data.index)),
                              index=data.index, name='impvolPut')

    for i in data.index:
        impvolput[i] = (calcvol(data.Expiration[i],
                                data.Strike[i],
                                today,
                                underlyingprice,
                                (data['Bid.1'][i] + data['Ask.1'][i]) / 2.0, Nag_Put))

    print 'Calculated Implied Vol for %i Puts' % len(data.index)

    data = data.join(impvolput)
    fig = plt.figure()
    fig.subplots_adjust(hspace=.4, wspace=.3)

    # Encode graph layout: 3 rows, 3 columns, 1 is first graph.
    num = 331
    max_xticks = 4

    for date in dates:
        # add each subplot to the figure
        plot_year, plot_month = date.split()
        plot_date = (int(plot_year) - 13) * 365 + cumulative_month[plot_month]
        plot_call = data[(data.impvolCall > 0) &
                       (data.impvolCall < 1) &
                       (data.Expiration == plot_date) &
                       (data['Last Sale'] > 0)]
        plot_put = data[(data.impvolPut > 0) &
                        (data.impvolPut < 1) &
                        (data.Expiration == plot_date) &
                        (data['Last Sale.1'] > 0)]

        myfig = fig.add_subplot(num)
        xloc = plt.MaxNLocator(max_xticks)
        myfig.xaxis.set_major_locator(xloc)
        myfig.set_title('Expiry: %s 20%s' % (plot_month, plot_year))
        myfig.plot(plot_call.Strike, plot_call.impvolCall, 'pr', label='call')
        myfig.plot(plot_put.Strike, plot_put.impvolPut, 'p', label='put')
        myfig.legend(loc=1, numpoints=1, prop={'size': 10})
        myfig.set_ylim([0,1])
        myfig.set_xlabel('Strike Price')
        myfig.set_ylabel('Implied Volatility')
        num += 1

    plt.suptitle('Implied Volatility for %s Current Price: %s Date: %s' %
                 (company, underlyingprice, qd_date))
    plt.show()
 
if __name__ == "__main__":
    main()
