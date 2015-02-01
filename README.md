# Implied Volatility Using Python's Pandas Package

This code is to accompany the corresponding Numerical Algorithms Group (NAG) blog post.

> http://blog.nag.com/2013/10/implied-volatility-using-pythons-pandas.html

This script uses options data downloaded from the [CBOE] in csv format. Be sure to download data during CBOE trading hours to ensure the data is not null. To run type 
'''sh
$ python implied_volatility.py QuoteData.dat
'''

This script has been tested with and requires the following packages/libraries:

  - Python 2.7 or 3.4
  - numpy 1.9.0
  - pandas 0.14.1
  - matplotlib 1.4.2
  - nag4py-23.0rc1

A [NAG C Library] license is required to run the script. To obtain a free 30-day trial license, email support@nag.co.uk

![Alt text](/pictures/VolatilityCurves.png?raw=true "Volatility Curves for AAPL")
![Alt text](/pictures/VolatilitySurface.png?raw=true "Volatility Surface for AAPL")

[NAG C Library]: http://www.nag.com/numeric/CL/CLdescription.asp
[CBOE]: http://www.cboe.com/delayedquote/QuoteTableDownload.aspx
