package net.imglib2.imglib2access;

import ij.ImageJ;
import ij.ImagePlus;

public class Code {

	static public ImgLib2Access filterBox_2D_NonSeparable( final ImgLib2Access in, final long length )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		final ImgLib2Access out = new ImgLib2Access( nx, ny );
		for ( long y = 0; y < ny; ++y )
			for ( long x = 0; x < nx; ++x )
			{
				final ImgLib2Access block = in.getNeighborhood( x, y, length, length );
				double sum = 0.0;
				for ( long i = 0; i < length; ++i )
					for ( long j = 0; j < length; ++j )
						sum += block.getPixel( i, j );
				
				out.putPixel( sum / ( length * length ), x, y );
			}
		return out;
	}

	static public ImgLib2Access filterBox_2D_Separable( ImgLib2Access in, final long length )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		ImgLib2Access out = new ImgLib2Access( nx, ny );
		for ( long y = 0; y < ny; ++y )
		{
			final ImgLib2Access rowin = in.getRow( y );
			final ImgLib2Access rowout = out.getRow( y );
			box1D( rowin, rowout, length );
		}
		in = out;
		out = new ImgLib2Access( nx, ny );
		for ( long x = 0; x < nx; ++x )
		{
			final ImgLib2Access colin = in.getColumn( x );
			final ImgLib2Access colout = out.getColumn( x );
			box1D( colin, colout, length );
		}
		return out;
	}

	static public ImgLib2Access filterBox_2D_Recursive( ImgLib2Access in, final long length )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		ImgLib2Access out = new ImgLib2Access( nx, ny );
		for ( long y = 0; y < ny; y++ )
		{
			final ImgLib2Access rowin = in.getRow( y );
			final ImgLib2Access rowout = out.getRow( y );
			box1Drecursive( rowin, rowout, length );
		}
		in = out;
		out = new ImgLib2Access( nx, ny );
		for ( long x = 0; x < nx; ++x )
		{
			final ImgLib2Access colin = in.getColumn( x );
			final ImgLib2Access colout = out.getColumn( x );
			box1Drecursive( colin, colout, length );
		}
		return out;
	}

	static protected void box1D( final ImgLib2Access u, final ImgLib2Access v, final long length )
	{
		final long n = u.getWidth();
		final long h = length / 2;
		for ( long k = 0; k < n; k++ )
		{
			v.putPixel( 0, k );
			for ( long i = -h; i <= h; i++ )
			{
				final long index = i + k;
				v.putPixel( v.getPixel( k ) + u.getPixel( index ), k );
			}
			v.putPixel( v.getPixel( k ) / length, k );
		}
	}

	static protected void box1Drecursive( final ImgLib2Access u, final ImgLib2Access v, final long length )
	{
		final long n = u.getWidth();
		final long h = length / 2;
		v.putPixel( 0, 0 );
		for ( long i = 0; i < length; ++i )
		{
			final long index = Math.abs( i - h );
			v.putPixel( v.getPixel( 0 ) + u.getPixel( index ), 0 );
		}
		v.putPixel( v.getPixel( 0 ) / length, 0 );
		long k1, k2;
		for ( long k = 1; k < n; ++k )
		{
			k1 = k - h - 1;
			if ( k1 < 0 ) k1 = -k1;
			k2 = k + h;
			if ( k2 >= n ) k2 = 2 * ( n - 1 ) - k2;
			v.putPixel( v.getPixel( k - 1 ) + ( u.getPixel( k2 ) - u.getPixel( k1 ) ) / length , k );
		}
	}

	static public ImgLib2Access filterA_2D( final ImgLib2Access in )
	{
		return filter3_2D( in, ImgLib2Access.row( -1.0, 0.0, 1.0 ), ImgLib2Access.row( 1.0, 1.0, 1.0 ) );
	}

	static public ImgLib2Access filterB_2D( final ImgLib2Access in )
	{
		return filter3_2D( in, ImgLib2Access.row( 0.25, 0.5, 0.25 ), ImgLib2Access.row( 0.25, 0.5, 0.25 ) );
	}

	static public ImgLib2Access filterC_2D( final ImgLib2Access in )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		final ImgLib2Access out = new ImgLib2Access( nx, ny );
		final ImgLib2Access coef = new ImgLib2Access( 5, 5 );
		coef.putPixels(
				-1.0, -1.0, 0.0, -1.0, -1.0,
				0.0, -1.0, 1.0, -1.0, 0.0,
				0.0, 0.0, 4.0, 0.0, 0.0,
				1.0, 1.0, 1.0, 1.0, 1.0,
				0.0, 1.0, 1.0, 1.0, 0.0 );
		for ( long x = 0; x < nx; x++ )
			for ( long y = 0; y < ny; y++ )
			{
				final ImgLib2Access block = in.getNeighborhood( x, y, 5, 5 );
				double sum = 0.0;
				for ( long k = 0; k < 5; k++ )
					for ( long l = 0; l < 5; l++ )
						sum += coef.getPixel( k, l ) * block.getPixel( k, l );
				out.putPixel( sum, x, y );
			}
		return out;
	}

	static public ImgLib2Access filterD_2D( final ImgLib2Access in )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		final long length = 5;
		final double a = -1.0 / ( length * length );
		final ImgLib2Access coef = new ImgLib2Access( length, length );
		coef.putPixels(
				a, a, a, a, a,
				a, a, a, a, a,
				a, a, 2.0 + a, a, a,
				a, a, a, a, a,
				a, a, a, a, a );
		final ImgLib2Access out = new ImgLib2Access( nx, ny );
		for ( long x = 0; x < nx; x++ )
			for ( long y = 0; y < ny; y++ )
			{
				final ImgLib2Access block = in.getNeighborhood( x, y, length, length );
				double sum = 0.0;
				for ( long k = 0; k < length; k++ )
					for ( long l = 0; l < length; l++ )
						sum += coef.getPixel( k, l ) * block.getPixel( k, l );
				out.putPixel( sum, x, y );
			}
		return out;
	}

	static protected ImgLib2Access filter3_2D( ImgLib2Access in, final ImgLib2Access maskRow, final ImgLib2Access maskCol )
	{
		final long nx = in.getWidth();
		final long ny = in.getHeight();
		ImgLib2Access out = new ImgLib2Access( nx, ny );
		for ( long y = 0; y < ny; y++ )
		{
			final ImgLib2Access rowin = in.getRow( y );
			final ImgLib2Access rowout = out.getRow( y );
			convolve3( rowin, rowout, maskRow );
		}
		in = out;
		out = new ImgLib2Access( nx, ny );
		for ( long x = 0; x < nx; x++ )
		{
			final ImgLib2Access colin = in.getColumn( x );
			final ImgLib2Access colout = out.getColumn( x );
			convolve3( colin, colout, maskCol );
		}
		return out;
	}

	static protected void convolve3( final ImgLib2Access vin, final ImgLib2Access vout, final ImgLib2Access mask )
	{
		assert vin.sizeEquals( vout ) : "Input and output sizes do not match.";
		
		final long n = vin.getWidth();
		vout.putPixel( mask.getPixel( 0 ) * vin.getPixel( 1 ) + mask.getPixel( 1 ) * vin.getPixel( 0 ) + mask.getPixel( 2 ) * vin.getPixel( 1 ), 0 );
		for ( long k = 1; k < n - 1; k++ )
			vout.putPixel( mask.getPixel( 0 ) * vin.getPixel( k - 1 ) + mask.getPixel( 1 ) * vin.getPixel( k ) + mask.getPixel( 2 ) * vin.getPixel( k + 1 ), k );
		
		vout.putPixel( mask.getPixel( 0 ) * vin.getPixel( n - 2 ) + mask.getPixel( 1 ) * vin.getPixel( n - 1 ) + mask.getPixel( 2 ) * vin.getPixel( n - 2 ), n - 1 );
	}

	static public ImgLib2Access localNormalization_2D( final ImgLib2Access in, final int length1, final int length2 )
	{
		final ImgLib2Access mean = filterBox_2D_Separable( in, length1 );
		final ImgLib2Access diff = in.duplicate();
		diff.subtract( diff, mean );
		ImgLib2Access var = diff.duplicate();
		var.pow( 2.0 );
		var = filterBox_2D_Separable( var, length2 );
		var.sqrt();
		diff.divide( diff, var );
		return diff;
	}
	
	protected void fillRows( final ImgLib2Access ima )
	{
		
	}
	
	public static void main( final String[] args )
	{
		new ImageJ();
		
		final ImagePlus imp = new ImagePlus( "http://fly.mpi-cbg.de/~saalfeld/saalfeld.jpg" );
		
		//final ImagePlus imp = IJ.createHyperStack( "3D", 200, 100, 3, 10, 5, 8 );
		//imp.putPixel( 255, 2, 4, 0 );
		
		imp.show();
		
		final ImgLib2Access ia = new ImgLib2Access( imp );
		filterBox_2D_NonSeparable( ia, 11 ).show( "box filter 2d non-separable" );
		filterBox_2D_Separable( ia, 11 ).show( "box filter 2d separable" );
		filterBox_2D_Recursive( ia, 11 ).show( "box filter 2d recursive" );
		filterA_2D( ia ).show( "filter A 2d" );
		filterB_2D( ia ).show( "filter B 2d" );
		filterC_2D( ia ).show( "filter C 2d" );
		filterD_2D( ia ).show( "filter C 2d" );
		localNormalization_2D( ia, 11, 81 ).show( "local normalized 2d" );
		
		//ib.createFloatImagePlus().show();
	}
}
