## Implied Volatility

An example script that will run on data from the [Chicago Board of Options Exchange](http://www.cboe.com/delayedquote/QuoteTableDownload.aspx) (CBOE).

The program uses the pandas package to easily store and manipulate the data via DataFrames.
In the script, you'll notice pandas takes care of the many data processing functions and
then calls the NAG Library for more complex analysis. The program will automatically read in
the options data, calculate implied volatility for the call and put options, and plot the
volatility curves and surface.

The above code can be run as follows (given that you have pandas, matplotlib, nag4py, and ctypes):

## Dependencies

 - NAG C Library
 - nag4py
 - pandas
 - matplotlib

If not already, the NAG C library also needs to be on your `LD_LIBRARY_PATH`.
