/**
 * License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package net.imglib2.imglib2access;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.converter.Converter;
import net.imglib2.display.projector.RandomAccessibleProjector2D;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.ByteArray;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgs;
import net.imglib2.img.planar.PlanarImgs;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.interpolation.randomaccess.LanczosInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.outofbounds.Bounded;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

/**
 * ImgLib2Access is a simplification layer on top of ImgLib2. It's aim is to
 * simplify teaching students in image processing beginners courses by removing
 * complexity. ImgLib2Access is a translation and extension of ImageAccess, a
 * similar interface to ImageJ by Sage and Unser (2001)
 *
 * BibTeX:
 * <pre>
 * &#64;inproceedings{958110,
 *   author    = {Sage, D. and Unser, M.},
 *   booktitle = {2001 International Conference on Image Processing, 2001},
 *   title     = {Easy Java programming for teaching image-processing},
 *   year      = {2001},
 *   volume    = {3},
 *   pages     = {298--301},
 *   doi       = {10.1109/ICIP.2001.958110},
 * }
 * </pre>
 *
 * The data are stored in an ImgLib2 {@link RandomAccessible} of
 * {@link RealType}s. {@link RandomAccessibleInterval}s are virtually extended
 * using one of ImgLib2's extension strategies and interpolated using one of
 * ImgLib2's interpolation strategies.
 *
 * @author Stephan Saalfeld <saalfelds@janelia.hhmi.org>
 */
public class ImgLib2Access
{
	final static protected class Infinite implements Bounded
	{
		public boolean isOutOfBounds()
		{
			return false;
		}
	}

	final static protected class FinalRealFloatConverter implements Converter< RealType< ? >, FloatType >
	{
		public void convert( final RealType< ? > input, final FloatType output )
		{
			output.set( input.getRealFloat() );
		}
	}

	final static protected class FinalRealUnsignedByteConverter implements Converter< RealType< ? >, UnsignedByteType >
	{
		final static protected int roundPositive( final double a )
		{
			return ( int )( a + 0.5 );
		}

		public void convert( final RealType< ? > input, final UnsignedByteType output )
		{
			output.set( Math.min( 255, roundPositive( Math.max( 0, input.getRealDouble() ) ) ) );
		}
	}

	/**
	 * Boundary extension methods
	 */
	static public enum Extension
	{
		ZERO, NAN, PERIODIC, MIRROR1, MIRROR2, RANDOM
	}

	/**
	 * Interpolation methods
	 */
	static public enum Interpolation
	{
		NN, NLINEAR, NLANCZOS
	}

	final protected RandomAccessibleInterval< RealType< ? > > srcInterval;
	final protected RandomAccessible< RealType< ? > > src;
	final protected RandomAccess< RealType< ? > > randomAccess;
	final protected Bounded bounded;
	final protected RealRandomAccess< RealType< ? > > realRandomAccess;


	/* Helpers for hiding generic type information */

	@SuppressWarnings( { "rawtypes", "unchecked" } )
	final static protected RandomAccessible< RealType< ? > > extend(
			final RandomAccessibleInterval< ? extends RealType< ? > > src, final Extension extension )
	{
		switch ( extension )
		{
		case PERIODIC:
			return ( RandomAccessible< RealType< ? > > ) Views.extendPeriodic( src );
		case MIRROR1:
			return ( RandomAccessible< RealType< ? > > ) Views.extendMirrorSingle( src );
		case MIRROR2:
			return ( RandomAccessible< RealType< ? > > ) Views.extendMirrorDouble( src );
		case RANDOM:
			return ( RandomAccessibleInterval ) extendRandom( src, 0, 1 );
		case NAN:
			{
				final RealType< ? > t = src.randomAccess().get().createVariable();
				t.setReal( Double.NaN );
				return ( RandomAccessible< RealType< ? > > ) Views.extendValue( ( RandomAccessibleInterval ) src, ( RealType ) t );
			}
		default:
			{
				final RealType< ? > t = src.randomAccess().get().createVariable();
				t.setZero();
				return ( RandomAccessible< RealType< ? > > ) Views.extendValue( ( RandomAccessibleInterval ) src, ( RealType ) t );
			}
		}
	}

	final static protected RandomAccessible< RealType< ? > > src(
			final RandomAccessible< RealType< ? > > src )
	{
		if ( RandomAccessibleInterval.class.isInstance( src ) ) return Views.extendMirrorSingle( ( RandomAccessibleInterval< RealType< ? > > ) src );
		else return src;
	}

	@SuppressWarnings( "unchecked" )
	final static protected < T extends RealType< T > > RandomAccessible< T > extendRandom(
			final RandomAccessibleInterval< ? extends RealType< ? > > src, final double min, final double max )
	{
		return Views.extendRandom( ( RandomAccessibleInterval< T > ) src, 0, 1 );
	}

	@SuppressWarnings( "unchecked" )
	final static protected < T extends RealType< T > > RealRandomAccess< T > interpolate(
			final RandomAccessible< ? extends RealType< ? > > src, final Interpolation interpolation )
	{
		final InterpolatorFactory< T, RandomAccessible< T > > factory;
		switch ( interpolation )
		{
		case NLINEAR:
			factory = new NLinearInterpolatorFactory< T >();
			break;
		case NLANCZOS:
			factory = new LanczosInterpolatorFactory< T >( 3, false );
			break;
		default:
			factory = new NearestNeighborInterpolatorFactory< T >();
		}

		final RealRandomAccessible< T > interpolant = Views.interpolate( ( RandomAccessible< T > )src, factory );
		return interpolant.realRandomAccess();
	}

	final static protected Bounded bounded( final RandomAccess< ? > randomAccess )
	{
		if ( Bounded.class.isInstance( randomAccess ) )
			return ( Bounded )randomAccess;
		else
			return new Infinite();
	}


	/* Constructors */

	@SuppressWarnings( "unchecked" )
	public ImgLib2Access( final RandomAccessible< RealType< ? > > src, final Interpolation interpolation, final long... size )
	{
		assert src.numDimensions() == size.length: "Numbers of dimensions do not match.";

		final FinalInterval interval = new FinalInterval( size );
		srcInterval = Views.interval( src, interval );
		this.src = src( src );
		randomAccess = src.randomAccess();
		bounded = bounded( randomAccess );
		realRandomAccess = ( RealRandomAccess< RealType< ? > > )( Object )interpolate( src, interpolation );
	}

	public ImgLib2Access( final RandomAccessible< RealType< ? > > src, final long... size )
	{
		this( src, Interpolation.NLINEAR, size );
	}

	@SuppressWarnings( "unchecked" )
	public ImgLib2Access( final RandomAccessibleInterval< RealType< ? > > srcInterval, final Extension extension, final Interpolation interpolation )
	{
		this.srcInterval = srcInterval;
		src = extend( srcInterval, extension );
		randomAccess = src.randomAccess();
		bounded = bounded( randomAccess );
		realRandomAccess = ( RealRandomAccess< RealType< ? > > )( Object )interpolate( src, interpolation );
	}

	public ImgLib2Access( final RandomAccessibleInterval< RealType< ? > > src, final Interpolation interpolation )
	{
		this( src, Extension.MIRROR1, interpolation );
	}

	public ImgLib2Access( final RandomAccessibleInterval< RealType< ? > > src, final Extension extension )
	{
		this( src, extension, Interpolation.NLINEAR );
	}

	public ImgLib2Access( final RandomAccessibleInterval< RealType< ? > > src )
	{
		this( src, Extension.MIRROR1, Interpolation.NLINEAR );
	}

	@SuppressWarnings( { "rawtypes", "unchecked" } )
	public ImgLib2Access( final ImagePlus imp, final Extension extension, final Interpolation interpolation ) throws RuntimeException
	{
		if ( imp.getType() != ImagePlus.COLOR_RGB )
		{
			final ImagePlusImg img = ( ImagePlusImg ) ImagePlusImgs.from( imp );
			srcInterval = img;
			src = extend( img, extension );
			randomAccess = src.randomAccess();
			bounded = bounded( randomAccess );
			realRandomAccess = ( RealRandomAccess< RealType< ? > > )( Object )interpolate( src, interpolation );
		}
		else throw new RuntimeException( "No int mapped ARGB color images supported, please make a composite." );

	}

	public ImgLib2Access( final ImagePlus imp, final Interpolation interpolation ) throws RuntimeException
	{
		this( imp, Extension.MIRROR1, interpolation );
	}

	public ImgLib2Access( final ImagePlus imp, final Extension extension ) throws RuntimeException
	{
		this( imp, extension, Interpolation.NLINEAR );
	}

	public ImgLib2Access( final ImagePlus imp ) throws RuntimeException
	{
		this( imp, Extension.MIRROR1, Interpolation.NLINEAR );
	}


	@SuppressWarnings( "unchecked" )
	public ImgLib2Access( final Extension extension, final Interpolation interpolation, final long... size )
	{
		/* accumulate */
		long n = size[ 0 ];
		for ( int i = 1; i < size.length; ++i )
			n *= size[ i ];

		final Img< DoubleType > img;
		if ( n <= 2147483648L ) img = ArrayImgs.doubles( size );
		else if ( size.length > 1 && size[ 0 ] * size[ 1 ] < 2147483648L ) img = PlanarImgs.doubles( size );
		else
		{
			/* suboptimal max cell-size, please fix if you really need it */
			final int[] cellSize = new int[ size.length ];
			final int maxSquareSize = ( int ) Math.pow( 2147483647, -size.length );
			for ( int i = 0; i < size.length; ++i )
				cellSize[ i ] = ( int ) Math.min( size[ i ], maxSquareSize );
			final CellImgFactory< DoubleType > factory = new CellImgFactory< DoubleType >( cellSize );
			img = factory.create( size, new DoubleType() );
		}

		srcInterval = ( RandomAccessibleInterval< RealType< ? > > ) ( Object ) img;
		src = extend( img, extension );
		randomAccess = src.randomAccess();
		bounded = bounded( randomAccess );
		realRandomAccess = ( RealRandomAccess< RealType< ? > > )( Object )interpolate( src, interpolation );
	}

	public ImgLib2Access( final Interpolation interpolation, final long... size )
	{
		this( Extension.MIRROR1, interpolation, size );
	}

	public ImgLib2Access( final Extension extension, final long... size )
	{
		this( extension, Interpolation.NLINEAR, size );
	}

	public ImgLib2Access( final long... size )
	{
		this( Extension.MIRROR1, Interpolation.NLINEAR, size );
	}

	static public ImgLib2Access row( final double... pixels )
	{
		final ImgLib2Access ila = new ImgLib2Access( pixels.length );
		ila.putPixels( pixels );
		return ila;
	}


	/* Methods */

	/**
	 * Return the number of dimensions of this {@link ImgLib2Access}.
	 *
	 * @return number of dimensions
	 */
	public int n()
	{
		return src.numDimensions();
	}

	/**
	 * Return the size of in some dimension of the {@link ImgLib2Access}.
	 *
	 * @return the image width
	 */
	public long[] getSize()
	{
		final long[] size = new long[ src.numDimensions() ];
		srcInterval.dimensions( size );
		return size;
	}

	/**
	 * Return the size of in some dimension of the {@link ImgLib2Access}.
	 *
	 * @return the image width
	 */
	public long getSize( final int d )
	{
		assert d >= 0 && d < src.numDimensions(): d + " is beyond the available dimensions.";

		return srcInterval.dimension( d );
	}

	/**
	 * Return the width of the image.
	 *
	 * @return the image width
	 */
	public long getWidth()
	{
		return srcInterval.dimension( 0 );
	}

	/**
	 * Return the height of the image (if it has one).
	 *
	 * @return the image height
	 */
	public long getHeight()
	{
		return getSize( 1 );
	}

	/**
	 * Returns true if the dimensions of this {@link ImgLib2Access} are equal
	 * to those of another {@link ImgLib2Access}.
	 *
	 * @param other
	 * 		the other {@link ImgLib2Access}
	 * @return true or false
	 */
	public boolean sizeEquals( final ImgLib2Access other )
	{
		boolean equals = true;
		for ( int d = 0; d < src.numDimensions(); ++d )
			equals &= srcInterval.dimension( d ) == other.srcInterval.dimension( d );
		return equals;
	}

	/**
	 * Return the maximum value of {@link ImgLib2Access}.
	 *
	 * @return the maximum value
	 */
	public double getMaximum()
	{
		final IterableInterval< RealType< ? > > iterable = Views.iterable( srcInterval );
		double max = -Double.MAX_VALUE;
		for ( final RealType< ? > t : iterable )
		{
			final double v = t.getRealDouble();
			if ( v > max ) max = v;
		}
		return max;
	}

	/**
	 * Return the minimum value of {@link ImgLib2Access}.
	 *
	 * @return the minimum value
	 */
	public double getMinimum()
	{
		final IterableInterval< RealType< ? > > iterable = Views.iterable( srcInterval );
		double min = Double.MAX_VALUE;
		for ( final RealType< ? > t : iterable )
		{
			final double v = t.getRealDouble();
			if ( v < min ) min = v;
		}
		return min;
	}

	/**
	 * Return the mean value of {@link ImgLib2Access}.
	 *
	 * @return the mean value
	 */
	public double getMean()
	{
		final IterableInterval< RealType< ? > > iterable = Views.iterable( srcInterval );
		double mean = 0.0;
		for ( final RealType< ? > t : iterable )
			mean += t.getRealDouble();
		mean /= ( double ) iterable.size();
		return mean;
	}

	final static protected ValuePair< ImageProcessor, RandomAccessibleProjector2D< RealType< ? >, ? > > createFloatProjector(
			final RandomAccessible< RealType< ? > > src, final int w, final int h )
	{
		final FloatProcessor fp = new FloatProcessor( w, h );
		final ArrayImg< FloatType, FloatArray > img = ArrayImgs.floats( ( float[] )fp.getPixels(), new long[]{ w, h } );
		final RandomAccessibleProjector2D< RealType< ? >, FloatType > projector =
		        new RandomAccessibleProjector2D< RealType< ? >, FloatType >( 0, 1, src, img, new FinalRealFloatConverter() );

		return new ValuePair< ImageProcessor, RandomAccessibleProjector2D< RealType< ? >, ? > >( fp, projector );
	}

	final static protected ValuePair< ImageProcessor, RandomAccessibleProjector2D< RealType< ? >, ? > > createUnsignedByteProjector(
			final RandomAccessible< RealType< ? > > src, final int w, final int h )
	{
		final ByteProcessor fp = new ByteProcessor( w, h );
		final ArrayImg< UnsignedByteType, ByteArray > img = ArrayImgs.unsignedBytes( ( byte[] )fp.getPixels(), new long[]{ w, h } );
		final RandomAccessibleProjector2D< RealType< ? >, UnsignedByteType > projector =
		        new RandomAccessibleProjector2D< RealType< ? >, UnsignedByteType >( 0, 1, src, img, new FinalRealUnsignedByteConverter() );

		return new ValuePair< ImageProcessor, RandomAccessibleProjector2D< RealType< ? >, ? > >( fp, projector );
	}


	/**
	 * Create an ImagePlus from the pixel data (works only for 1 &lt;
	 * <em>n</em> &lt; 6). The double values of the pixel are simply casted to
	 * float.
	 *
	 * @return an ImagePlus
	 */
	public ImagePlus createImagePlus( final int type )
	{
		assert
			src.numDimensions() > 1 &&
			src.numDimensions() < 6 &&
			getWidth() * getHeight() < 2147483648L: "ImagePlus supports only 2 to 5 dimensions with less than 2^31 pixels.";

		final int n = src.numDimensions();
		final int w = ( int )getWidth();
		final int h = ( int )getHeight();
		final ImageStack stack = new ImageStack( w, h );

		final ValuePair< ImageProcessor, RandomAccessibleProjector2D< RealType< ? >, ? > > pair;
		switch ( type )
		{
		case ImagePlus.GRAY8:
			pair = createUnsignedByteProjector( src, w, h );
			break;
		default:
			pair = createFloatProjector( src, w, h );
		}

		int d = 2;
		do
		{
			pair.getB().map();
			stack.addSlice( pair.a.duplicate() );

			for ( d = 2; d < n; ++d )
			{
				pair.b.fwd( d );
				if ( pair.b.getLongPosition( d ) < getSize( d ) )
					break;
				else
					pair.b.setPosition( 0, d );
			}
		}
		while ( d < n && pair.b.getLongPosition( n - 1 ) < getSize( n - 1 ) );

		final ImagePlus imp = new ImagePlus( "", stack );
		if ( n > 2 )
		{
			imp.setOpenAsHyperStack( true );
			final int c = ( int )getSize( 2 ), s, f;
			if ( n > 3 )
			{
				s = ( int )getSize( 3 );
				if ( n > 4 )
					f = ( int )getSize( 4 );
				else
					f = 1;
			}
			else
			{
				s = 1;
				f = 1;
			}
			imp.setDimensions( c, s, f );
		}
		return imp;
	}

	/**
	 * Create a FloatImagePlus from the pixel data (works only for 1 &lt;
	 * <em>n</em> &lt; 6). The double values of the pixel are simply casted to
	 * float.
	 *
	 * @return the FloatProcessor
	 */
	public ImagePlus createFloatImagePlus()
	{
		return createImagePlus( ImagePlus.GRAY32 );
	}

	/**
	 * Create a ByteImagePlus from the pixel data (works only for 1 &lt;
	 * <em>n</em> &lt; 6). The double values of the pixel are casted and cropped to
	 * [0-255].
	 *
	 * @return the FloatProcessor
	 */
	public ImagePlus createByteImagePlus()
	{
		return createImagePlus( ImagePlus.GRAY8 );
	}

	/**
	 * Create a new {@link ImgLib2Access} object by deep copying the current
	 * {@link ImgLib2Access} object's pixels.
	 *
	 * @return a new {@link ImgLib2Access} object
	 **/
	public ImgLib2Access duplicate()
	{
		final ImgLib2Access ia = new ImgLib2Access( getSize() );

		/* safe copy by flat iteration order */
		final Cursor< RealType< ? > > c1 = Views.flatIterable( srcInterval ).cursor();
		final Cursor< RealType< ? > > c2 = Views.flatIterable( ia.srcInterval ).cursor();

		while ( c1.hasNext() )
			c2.next().setReal( c1.next().getRealDouble() );

		return ia;
	}


	/**
	 * An ImageAccess object calls this method for getting the gray level of a
	 * selected pixel.
	 *
	 * @param x
	 *            input, the integer x-coordinate of a pixel
	 * @param y
	 *            input, the integer y-coordinate of a pixel
	 * @return the gray level of the pixel (double)
	 */
	public double getPixel( final long... position )
	{
		assert src.numDimensions() <= position.length : "Numbers of dimensions do not match.";

		randomAccess.setPosition( position );
		return randomAccess.get().getRealDouble();
	}

	/**
	 * An ImageAccess object calls this method for getting the gray level of a
	 * selected pixel at an interpolated real location.
	 *
	 * @param x
	 *            input, the double x-coordinate of a pixel
	 * @param y
	 *            input, the double y-coordinate of a pixel
	 * @return the gray level of the pixel (double)
	 */
	public double getInterpolatedPixel( final double... position )
	{
		assert src.numDimensions() <= position.length : "Numbers of dimensions do not match.";

		realRandomAccess.setPosition( position );
		return realRandomAccess.get().getRealDouble();
	}


	final static protected boolean columnInRange( final Interval interval, final long... position )
	{
		assert position.length >= interval.numDimensions() - 1 : "Numbers of dimensions do not match.";

		boolean is = position[ 0 ] >= interval.min( 0 ) && position[ 0 ] <= interval.max( 0 );
		for ( int d = 2; d < interval.numDimensions() && is; ++d )
		{
			final long p = position[ d - 1 ];
			is &= p >= interval.min( d ) && p <= interval.max( d );
		}
		return is;
	}

	/**
	 * An {@link ImgLib2Access} object calls this method for getting a whole
	 * column of the image.
	 *
	 * @param x[, z, t, ...] the coordinate of the column, skipping the second dimension
	 *
	 * @return	the column as an ImgLib2Access
	 */
	public ImgLib2Access getColumn( final long... position )
	{
		assert columnInRange( srcInterval, position ) : "Position out of range.";

		RandomAccessibleInterval< RealType< ? > > column = Views.hyperSlice( srcInterval, 0, position[ 0 ] );
		for ( int d = 2; d < src.numDimensions(); ++d )
			column = Views.hyperSlice( column, 0, position[ d - 1 ] );

		return new ImgLib2Access( column );
	}


	final static protected boolean rowInRange( final Interval interval, final long... position )
	{
		assert position.length >= interval.numDimensions() - 1 : "Numbers of dimensions do not match.";

		boolean is = position[ 0 ] >= interval.min( 1 ) && position[ 0 ] <= interval.max( 1 );
		for ( int d = 2; d < interval.numDimensions() && is; ++d )
		{
			final long p = position[ d - 1 ];
			is &= p >= interval.min( d ) && p <= interval.max( d );
		}
		return is;
	}

	/**
	 * An {@link ImgLib2Access} object calls this method for getting a whole
	 * row of the image.
	 *
	 * @param y[, z, t, ...] the coordinate of the row, skipping the first dimension
	 *
	 * @return	the row as an ImgLib2Access
	 */
	public ImgLib2Access getRow( final long... position )
	{
		assert rowInRange( srcInterval, position ) : "Position out of range.";

		RandomAccessibleInterval< RealType< ? > > row = Views.hyperSlice( srcInterval, 1, position[ 0 ] );
		for ( int d = 2; d < src.numDimensions(); ++d )
			row = Views.hyperSlice( row, 0, position[ d - 1 ] );

		return new ImgLib2Access( row );
	}


	final static protected boolean kRowInRange( final Interval interval, final int k, final long... position )
	{
		assert
			position.length >= interval.numDimensions() - 1 &&
			k >= 0 &&
			k < interval.numDimensions() : "Numbers of dimensions do not match.";

		boolean is = true;
		for ( int d = 0; d < k && is; ++d )
		{
			is &= position[ d ] >= interval.min( d ) && position[ d ] <= interval.max( d );
		}
		for ( int d = k + 1; d < interval.numDimensions() && is; ++d )
		{
			final long p = position[ d - 1 ];
			is &= p >= interval.min( d ) && p <= interval.max( d );
		}
		return is;
	}

	/**
	 * Get a whole pixel row in dimension <em>k</em> at a specified position.
	 *
	 * @param <em>x</em>[, <em>y</em>, ...] the coordinates of the
	 * 		<em>k</em>-row, skipping the <em>k</em>-th dimension
	 *
	 * @return	the <em>k</em>-row as an {@link ImgLib2Access}
	 */
	public ImgLib2Access getKRow( final int k, final long... position )
	{
		assert kRowInRange( srcInterval, k, position ) : "Position out of range.";

		RandomAccessibleInterval< RealType< ? > > kRow = srcInterval;

		/* before k */
		for ( int d = 0; d < k; ++d )
			kRow = Views.hyperSlice( kRow, 0, position[ d ] );

		/* after k */
		for ( int d = k + 1; d < src.numDimensions(); ++d )
			kRow = Views.hyperSlice( kRow, 1, position[ d - 1 ] );

		return new ImgLib2Access( kRow );
	}

	/**
	 * Get a neighborhood around a pixel position.
	 *
	 * @param parameters
	 * 			the integer x[, y, ...]-coordinates of a selected central
	 *			pixel and the width[, height, ...]-dimension of the
	 *			neighborhood.
	 * @return
	 *			neighborhood
	 */
	public ImgLib2Access getNeighborhood( final long... parameters )
	{
		final int n = src.numDimensions();

		assert parameters.length == 2 * n : "Dimensions do not match.";

		final long[] min = new long[ n ];
		final long[] max = new long[ n ];

		for ( int d = 0; d < n; ++d )
		{
			final long wd = parameters[ d + n ];;
			min[ d ] = parameters[ d ] - ( wd - 1 ) / 2;
			max[ d ] = min[ d ] + wd - 1;
		}

		return new ImgLib2Access( Views.offsetInterval( src, min, max ) );
	}

	/**
	 * Get a sub-image with the upper left corner in the coordinate
	 * (x[, y, ...]) with a size of (width[, height, ...]).
	 *
	 * @param parameters
	 *			the integer x[, y, ...]-coordinates of the sub-image and its
	 *			width[, height, ...]-dimension.
	 * @return
	 * 			sub-image
	 */
	public ImgLib2Access getSubImage( final long... parameters )
	{
		final int n = src.numDimensions();

		assert parameters.length == 2 * n : "Dimensions do not match.";

		final long[] min = new long[ n ];
		final long[] max = new long[ n ];

		for ( int d = 0; d < n; ++d )
		{
			final long wd = parameters[ d + n ];;
			min[ d ] = parameters[ d ];
			max[ d ] = min[ d ] + wd - 1;
		}

		return new ImgLib2Access( Views.offsetInterval( src, min, max ) );
	}

	/**
	 * Set a given gray value to a pixel at the given coordinate.  If the
	 * coordinate is outside the image boundaries, the set operation is
	 * ignored.
	 *
	 * @param x
	 *            input, the integer x-coordinate of a pixel
	 * @param y
	 *            input, the integer y-coordinate of a pixel
	 * @param value
	 *            input, a value of the gray level of the type double
	 */
	public void putPixel( final double value, final long... position )
	{
		assert src.numDimensions() <= position.length : "Numbers of dimensions do not match.";

		randomAccess.setPosition( position );
		if ( bounded.isOutOfBounds() )
			return;
		else
			randomAccess.get().setReal( value );
	}

	/**
	 * Copy a given {@link ImgLib2Access} into the current
	 * {@link ImgLib2Access} at a given offset.
	 *
	 * @param x
	 *            input, the integer x-coordinate of a pixel
	 * @param y
	 *            input, the integer y-coordinate of a pixel
	 * @param value
	 *            input, a value of the gray level of the type double
	 */
	public void putImage( final ImgLib2Access ima, final long... position )
	{
		assert
			ima.src.numDimensions() == src.numDimensions() &&
			src.numDimensions() <= position.length : "Numbers of dimensions do not match.";

		final long[] min = new long[ src.numDimensions() ];
		final long[] mix = new long[ min.length ];
		final long[] man = new long[ min.length ];
		final long[] max = new long[ min.length ];

		for ( int d = 0; d < min.length; ++d )
		{
			min[ d ] = Math.max( position[ d ], 0 );
			mix[ d ] = Math.min( position[ d ] + ima.srcInterval.dimension( d ) - 1, srcInterval.max( d ) );
			if ( mix[ d ] < min[ d ] )
				return;
			if ( position[ d ] < 0 )
				man[ d ] = -position[ d ];
			max[ d ] = mix[ d ] - min[ d ];
		}

		final RandomAccessibleInterval< RealType< ? > > a = Views.interval( ima.srcInterval, min, mix );
		final RandomAccessibleInterval< RealType< ? > > b = Views.interval( srcInterval, man, max );
		final Cursor< RealType< ? > > cA = Views.flatIterable( a ).cursor();
		final Cursor< RealType< ? > > cB = Views.flatIterable( b ).cursor();

		while ( cA.hasNext() )
			cB.next().setReal( cA.next().getRealDouble() );
	}

	/**
	 * Set the values of all pixels of the {@link ImgLib2Access}.
	 *
	 * @param pixels
	 *            pixel values
	 */
	public void putPixels( final double... pixels )
	{
		assert Views.iterable( srcInterval ).size() >= pixels.length : "Too many values.";

		int i = -1;
		for ( final RealType< ? > t : Views.flatIterable( srcInterval ) )
			t.setReal( pixels[ ++i ] );
	}

	/**
	 * Set a constant value to all pixels of the {@link ImgLib2Access}.
	 *
	 * @param constant
	 *            a constant value
	 */
	public void setConstant( final double constant )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( constant );
	}


	/**
	 * Stretches the contrast inside an image so that the gray levels are in the
	 * range 0 to 255.
	 *
	 * @return this
	 */
	public ImgLib2Access normalizeContrast()
	{
		final double maxGoal = 255.0;

		// Search the min and max
		double minImage = getMinimum();
		final double maxImage = getMaximum();

		// Compute the parameter to rescale the gray levels
		double a;

		if ( minImage - maxImage == 0 )
		{
			a = 1.0;
			minImage = maxGoal / 2.0;
		}
		else
			a = maxGoal / ( maxImage - minImage );

		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( a * pixel.getRealDouble() - minImage );

		return this;
	}

	/**
	 * Display an image at a specific position (x, y).
	 *
	 * @param title of the window
	 * @param loc location
	 */
	public void show( final String title, final java.awt.Point loc )
	{
		final ImagePlus impResult = createFloatImagePlus();
		impResult.setTitle( title );
		impResult.getProcessor().resetMinAndMax();
		impResult.show();
		final ij.gui.ImageWindow window = impResult.getWindow();
		window.setLocation( loc.x, loc.y );
		impResult.show();
	}

	/**
	 * Display an image.
	 *
	 * @param title of the window
	 */
	public void show( final String title )
	{
		final ImagePlus impResult = createFloatImagePlus();
		impResult.setTitle( title );
		impResult.getProcessor().resetMinAndMax();
		impResult.show();
	}

	/**
	 * Display an image at a specific position (x, y) as unsigned gray values.
	 *
	 * @param title of the window
	 * @param loc location
	 */
	public void showAsBytes( final String title, final java.awt.Point loc )
	{
		final ImagePlus impResult = createByteImagePlus();
		impResult.setTitle( title );
		impResult.show();
		final ij.gui.ImageWindow window = impResult.getWindow();
		window.setLocation( loc.x, loc.y );
		impResult.show();
	}

	/**
	 * Display an image as a unsigned gray values.
	 *
	 * @param title of the window
	 */
	public void showAsBytes( final String title )
	{
		final ImagePlus impResult = createByteImagePlus();
		impResult.setTitle( title );
		impResult.show();
	}

	/**
	 * Compute the absolute values of an {@link ImgLib2Access}.
	 *
	 * @return this
	 */
	public ImgLib2Access abs()
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( Math.abs( pixel.getRealDouble() ) );

		return this;
	}

	/**
	 * Compute the square roots of an {@link ImgLib2Access}.
	 *
	 * @return this
	 */
	public ImgLib2Access sqrt()
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( Math.sqrt( pixel.getRealDouble() ) );

		return this;
	}

	/**
	 * Raise all pixels to the power a.
	 *
	 * @param a
	 *            input
	 *
	 * @return this
	 */
	public ImgLib2Access pow( final double a )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( Math.pow( pixel.getRealDouble(), a ) );

		return this;
	}

	/**
	 * Add a constant to each pixel.
	 *
	 * @param constant
	 *            a constant to be added
	 *
	 * @return this
	 */
	public ImgLib2Access add( final double constant )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( pixel.getRealDouble() + constant );

		return this;
	}

	/**
	 * Multiply each pixel with a constant.
	 *
	 * @param constant
	 *            a constant to be multiplied
	 *
	 * @return this
	 */
	public ImgLib2Access multiply( final double constant )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( pixel.getRealDouble() * constant );

		return this;
	}

	/**
	 * Subtract a constant from each pixel.
	 *
	 * @param constant
	 *            a constant to be subtracted
	 *
	 * @return this
	 */
	public ImgLib2Access subtract( final double constant )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( pixel.getRealDouble() - constant );

		return this;
	}

	/**
	 * Divide each pixel by a constant.
	 *
	 * @param constant
	 *            a constant to be subtracted
	 *
	 * @return this
	 */
	public ImgLib2Access divide( final double constant )
	{
		for ( final RealType< ? > pixel : Views.iterable( srcInterval ) )
			pixel.setReal( pixel.getRealDouble() - constant );

		return this;
	}

	/**
	 * Add two {@link ImgLib2Access}es into this {@link ImgLib2Access}.
	 *
	 * [this = im1 + im2]
	 *
	 * The resulting {@link ImgLib2Access} and the two operands should have the same size.
	 *
	 * @param im1
	 *            an {@link ImgLib2Access} to be added
	 * @param im2
	 *            an {@link ImgLib2Access} to be added
	 *
	 * @return this
	 */
	public ImgLib2Access add( final ImgLib2Access im1, final ImgLib2Access im2 )
	{
		assert ( sizeEquals( im1 ) && sizeEquals( im2 ) ) : "Sizes do not match.";

		final Cursor< RealType< ? > > a = Views.flatIterable( im1.srcInterval ).cursor();
		final Cursor< RealType< ? > > b = Views.flatIterable( im2.srcInterval ).cursor();
		final Cursor< RealType< ? > > c = Views.flatIterable( srcInterval ).cursor();

		while ( a.hasNext() )
			c.next().setReal( a.next().getRealDouble() + b.next().getRealDouble() );

		return this;
	}

	/**
	 * Multiply two {@link ImgLib2Access}es into this {@link ImgLib2Access}.
	 *
	 * [this = im1 * im2]
	 *
	 * The resulting {@link ImgLib2Access} and the two operands should have the same size.
	 *
	 * @param im1
	 *            an {@link ImgLib2Access} to be multiplied
	 * @param im2
	 *            an {@link ImgLib2Access} to be multiplied
	 *
	 * @return this
	 */
	public ImgLib2Access multiply( final ImgLib2Access im1, final ImgLib2Access im2 )
	{
		assert ( sizeEquals( im1 ) && sizeEquals( im2 ) ) : "Sizes do not match.";

		final Cursor< RealType< ? > > a = Views.flatIterable( im1.srcInterval ).cursor();
		final Cursor< RealType< ? > > b = Views.flatIterable( im2.srcInterval ).cursor();
		final Cursor< RealType< ? > > c = Views.flatIterable( srcInterval ).cursor();

		while ( a.hasNext() )
			c.next().setReal( a.next().getRealDouble() * b.next().getRealDouble() );

		return this;
	}

	/**
	 * Subtract two {@link ImgLib2Access}es into this {@link ImgLib2Access}.
	 *
	 * [this = im1 - im2]
	 *
	 * The resulting {@link ImgLib2Access} and the two operands should have the same size.
	 *
	 * @param im1
	 *            an {@link ImgLib2Access} to be subtracted
	 * @param im2
	 *            an {@link ImgLib2Access} to be subtracted
	 *
	 * @return this
	 */
	public ImgLib2Access subtract( final ImgLib2Access im1, final ImgLib2Access im2 )
	{
		assert ( sizeEquals( im1 ) && sizeEquals( im2 ) ) : "Sizes do not match.";

		final Cursor< RealType< ? > > a = Views.flatIterable( im1.srcInterval ).cursor();
		final Cursor< RealType< ? > > b = Views.flatIterable( im2.srcInterval ).cursor();
		final Cursor< RealType< ? > > c = Views.flatIterable( srcInterval ).cursor();

		while ( a.hasNext() )
			c.next().setReal( a.next().getRealDouble() - b.next().getRealDouble() );

		return this;
	}


	/**
	 * Divide two {@link ImgLib2Access}es into this {@link ImgLib2Access}.
	 *
	 * [this = im1 / im2]
	 *
	 * The resulting {@link ImgLib2Access} and the two operands should have the same size.
	 *
	 * @param im1
	 *            numerator
	 * @param im2
	 *            denominator
	 *
	 * @return this
	 */
	public ImgLib2Access divide( final ImgLib2Access im1, final ImgLib2Access im2 )
	{
		assert ( sizeEquals( im1 ) && sizeEquals( im2 ) ) : "Sizes do not match.";

		final Cursor< RealType< ? > > a = Views.flatIterable( im1.srcInterval ).cursor();
		final Cursor< RealType< ? > > b = Views.flatIterable( im2.srcInterval ).cursor();
		final Cursor< RealType< ? > > c = Views.flatIterable( srcInterval ).cursor();

		while ( a.hasNext() )
			c.next().setReal( a.next().getRealDouble() / b.next().getRealDouble() );

		return this;
	}
}
