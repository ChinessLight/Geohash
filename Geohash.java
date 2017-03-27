public class Geohash {

	public static final int MAX_PRECISION = 30;
	
	private static final long maxLngPrecision = 31;
	private static final long minHashPrecision31 = 1L << maxLngPrecision;

	private static final long maxHashPrecision30 = minHashPrecision31 - 1;
	private static final long maxHashPrecision31 = minHashPrecision31 | maxHashPrecision30;

	private static final double lngf = minHashPrecision31 / 360.0;
	private static final double latf = minHashPrecision31 / 180.0;
	
	public static final double SMALLEST_SPAN_DEGREES = 360.0 / minHashPrecision31;
	public static final double errorTol = SMALLEST_SPAN_DEGREES / 1000.0;
	
	private static long latHash(double lat, int precision) {
		if (precision < 0 || precision > MAX_PRECISION) {
			throw new IllegalArgumentException("invalid precision: " + precision);
		}
		long f = (long)((lat + 90.0) * latf);
		if (f > maxHashPrecision30) {
			return maxHashPrecision31 >>> (31 - precision);
		} else if (f <= 0L) {
			return minHashPrecision31 >>> (31 - precision);
		} else {
			return (f | minHashPrecision31) >>> (31 - precision);
		}
	}
	
	/**
	 * returns a longitude hash that is valid for the given latitude hash
	 */
	static long lngHash(long latHash) {
		if (latHash == 1L) { //special case
			return 1L;
		}
		long latHashMax = max(latHash);
		return lngHashFromMin(latHash, latHashMax - (latHashMax >> 1));
	}
	
	private static long lngHashFromMin(long latHash, long m) {
		return ((latHash ^ (max(latHash & (m >> 1)) | m)) << 3) | 4;
	}
	
	private static long lngHash(double lng, long latHash) {
		long aLngHash = lngHash(latHash);
		long f = (long) ((lng + 180.0) * lngf);
		return (f & maxHashPrecision30 | minHashPrecision31) >> invp(aLngHash);
	}
	
	static long valueOf(long latHash, double lng) {
		return (widen(latHash) << 1) | widen(lngHash(lng, latHash));
	}
	
	
	/**
	 * Returns the geohash which contains the specified coordinate at the specified precision. 
	 * Precision 0 represents the bounds of the entire earth and precision 30 represents a span 
	 * equal to approximately 2 centimeters. 
	 * <p>
	 * <strong>Note:</strong> this geohashing algorithm differentiates itself from
	 * the conventional geohashing algorithm due to the fact that the ratio of the latitudinal to longitudinal 
	 * distance (in meters, not degrees) across the center of a geohash is guaranteed to be greater than 
	 * or equal to 0.5, and less than 1.5. This makes it feasible to query regions using geohashes.
	 * @throws IllegalArgumentException if the precision is less than 0 or greater than 30.
	 */
	public static long valueOf(double lat, double lng, int pre) {
		return valueOf(latHash(lat, pre), lng);
	}
	
	/**
	 * Returns true if the geohash is valid.
	 */
	public static boolean isValid(long geohash) {
		long latHash = unwiden(geohash >> 1);
		if (latHash < 1 || latHash > maxHashPrecision30) {
			return false;
		}
		return max(unwiden(geohash)) == max(lngHash(latHash));
	}
	
	public static int getPrecision(long geohash) {
		return 31 - invp(unwiden(geohash >> 1));
	}
	
	public static double southernLatitude(long geohash) {
		long u = unwiden(geohash >> 1);
		return ((u << invp(u)) & maxHashPrecision30) / latf - 90.0;
	}
	
	public static double latitudeSpan(long geohash) {
		long u = unwiden(geohash >> 1);
		return 180.0 / (minHashPrecision31 >> invp(u));
	}
	
	public static double westernLongitude(long geohash) {
		long u = unwiden(geohash);
		return ((u << invp(u)) & maxHashPrecision30) / lngf - 180.0;
	}
	
	public static double longitudeSpan(long geohash) {
		long u = unwiden(geohash);
		return 360.0 / (minHashPrecision31 >> invp(u));
	}
	
	/**
	 * Returns true if the specified geohash intersects the region specified by the specified southwest and northeast coordinates.
	 * Results are undefined if the eastern or western longitudes of the given region are out of the range [-180, 180].
	 * @throws IllegalArgumentException if the southern latitude is greater than the northern latitude
	 */
	public static boolean intersects(long geohash, double southLat, double westLng, double northLat, double eastLng) {
		long latHash = unwiden(geohash >> 1);
		long lngHash = unwiden(geohash);
		int invpLat = invp(latHash);
		int invpLng = invp(lngHash);
		double slat = ((latHash << invpLat) & maxHashPrecision30) / latf - 90.0;
		double wlng = ((lngHash << invpLng) & maxHashPrecision30) / lngf - 180.0;
		double nlat = slat + 180.0 / (minHashPrecision31 >> invpLat);
		double elng = wlng + 360.0 / (minHashPrecision31 >> invpLng);
		
		if (northLat >= southLat) {
			if (southLat > nlat || northLat < slat) {
				return false;
			} 
		} else {
			throw new IllegalArgumentException(southLat + ">" + northLat);
		}
		
		if (eastLng >= westLng) {
			if (westLng == -180 && elng == 180 || wlng == -180 && eastLng == 180) {
				return true;
			}
			if (westLng > elng || eastLng < wlng) {
				return false;
			}
		} else {
			if (westLng > elng && eastLng < wlng) {
				return false;
			}
		}
		
		return true;
	}
	
	
	public static String toDebugString(long geohash) {
		double slat = southernLatitude(geohash);
		double wlng = westernLongitude(geohash);
		double nlat = slat + latitudeSpan(geohash);
		double elng = wlng + longitudeSpan(geohash);
		return "0b" + Long.toBinaryString(geohash) + 
				"<pre:" + getPrecision(geohash) + " SW:(" + (float)slat + "," + (float)wlng + 
				") NE:(" + (float)nlat + "," + (float)elng + ")>";
	}
	
	
	/**
	 * Shifts a geohash one spot to the north. If there are two geohashes in the northern spot, 
	 * returns the western one. If there are no geohashes in the northern spot, returns 0.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftNorth0(final long hash) {

		long oldLatHash = unwiden(hash >> 1);
		long latHashMax = max(oldLatHash);
		if (oldLatHash == latHashMax) { //this geohash doesn't have a north; we're at the north pole!
			return 0;
		}
		long latHashMin = latHashMax - (latHashMax >> 1);
		
		long newLatHash = oldLatHash + 1;
		
		long newLngHashMax = max(lngHashFromMin(newLatHash, latHashMin));
		long oldLngHashMax = max(lngHashFromMin(oldLatHash, latHashMin));
		
		long widenedLatHash = (widen(newLatHash) << 1);
		long widenedLngHash = (hash & 0x5555555555555555L);
				
		if (newLngHashMax == oldLngHashMax) {
			return widenedLatHash | widenedLngHash;
		} else if (newLngHashMax == (oldLngHashMax & newLngHashMax)) {
			return widenedLatHash | (widenedLngHash >> 2);
		} else {
			return widenedLatHash | (widenedLngHash << 2);
		}
	}
	
	
	/**
	 * Shifts a geohash one spot to the north. If there are two geohashes in the northern spot, 
	 * returns the eastern one. If there are no geohashes in the northern spot, returns 0.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftNorth1(final long hash) {

		long oldLatHash = unwiden(hash >> 1);
		long latHashMax = max(oldLatHash);
		if (oldLatHash == latHashMax) { //this geohash doesn't have a north; we're at the north pole!
			return 0;
		}
		long latHashMin = latHashMax - (latHashMax >> 1);

		long newLatHash = oldLatHash + 1;
		
		long newLngHashMax = max(lngHashFromMin(newLatHash, latHashMin));
		long oldLngHashMax = max(lngHashFromMin(oldLatHash, latHashMin));
		
		long widenedLatHash = (widen(newLatHash) << 1);
		long widenedLngHash = (hash & 0x5555555555555555L);
				
		if (newLngHashMax == oldLngHashMax) {
			return widenedLatHash | widenedLngHash;
		} else if (newLngHashMax == (oldLngHashMax & newLngHashMax)) {
			return widenedLatHash | (widenedLngHash >> 2);
		} else {
			return widenedLatHash | (widenedLngHash << 2) | 1;
		}
	}

	/**
	 * Shifts a geohash one spot to the south. If there are two geohashes in the southern spot, 
	 * returns the western one. If there are no geohashes in the southern spot, returns 0.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftSouth0(final long hash) {
		long oldLatHash = unwiden(hash >> 1);
		long latHashMax = max(oldLatHash);
		long latHashMin = latHashMax - (latHashMax >> 1);
		if (oldLatHash == latHashMin) { //this geohash doesn't have a south; we're at the south pole!
			return 0;
		}
		
		long newLatHash = oldLatHash - 1;
		
		long newLngHashMax = max(lngHashFromMin(newLatHash, latHashMin));
		long oldLngHashMax = max(lngHashFromMin(oldLatHash, latHashMin));
		
		long widenedLatHash = (widen(newLatHash) << 1);
		long widenedLngHash = (hash & 0x5555555555555555L);
				
		if (newLngHashMax == oldLngHashMax) {
			return widenedLatHash | widenedLngHash;
		} else if (newLngHashMax == (oldLngHashMax & newLngHashMax)) {
			return widenedLatHash | (widenedLngHash >> 2);
		} else {
			return widenedLatHash | (widenedLngHash << 2);
		}
	}
	
	/**
	 * Shifts a geohash one spot to the south. If there are two geohashes in the southern spot, 
	 * returns the eastern one. If there are no geohashes in the southern spot, returns 0.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftSouth1(final long hash) {
		long oldLatHash = unwiden(hash >> 1);
		long latHashMax = max(oldLatHash);
		long latHashMin = latHashMax - (latHashMax >> 1);
		if (oldLatHash == latHashMin) { //this geohash doesn't have a south; we're at the south pole!
			return 0;
		}
		
		long newLatHash = oldLatHash - 1;
		
		long newLngHashMax = max(lngHashFromMin(newLatHash, latHashMin));
		long oldLngHashMax = max(lngHashFromMin(oldLatHash, latHashMin));
		
		long widenedLatHash = (widen(newLatHash) << 1);
		long widenedLngHash = (hash & 0x5555555555555555L);
				
		if (newLngHashMax == oldLngHashMax) {
			return widenedLatHash | widenedLngHash;
		} else if (newLngHashMax == (oldLngHashMax & newLngHashMax)) {
			return widenedLatHash | (widenedLngHash >> 2);
		} else {
			return widenedLatHash | (widenedLngHash << 2) | 1;
		}
	}
	
	/**
	 * Returns the geohash that contains this geohash, or 0 if no geohash contains this geohash. 
	 * Results are undefined for invalid input geohashes.
	 */
	public static long zoomOut(final long hash) {
		if ((hash | 0b11111) == 0b11111) {
			return hash == 0b11 ? 0 : 0b11;
		} else {
			long widenedLatHash = hash & 0xaaaaaaaaaaaaaaaaL;
			long widenedLatHashMax = max(widenedLatHash);
			widenedLatHashMax |= (widenedLatHashMax >> 32);
			
			if (widenedLatHash == (widenedLatHashMax & 0xaaaaaaaaaaaaaaaaL) || widenedLatHash == (widenedLatHashMax - (widenedLatHashMax >> 1))) {
				return (widenedLatHash >> 2) | (hash & 0x5555555555555555L);
			} else {
				return hash >> 2;
			}
		}
	}
	
	private static long zoomInNorth1(final long hash) {
		if ((hash & 0x2000000000000000L) != 0) {
			return 0;
		}			
		long widenedLatHash = ((hash & 0xaaaaaaaaaaaaaaaaL) << 2) | 0b10;
		long widenedLatHashMax = max(widenedLatHash);
		widenedLatHashMax |= (widenedLatHashMax >> 32);
		
		if (widenedLatHash == (widenedLatHashMax & 0xaaaaaaaaaaaaaaaaL) || widenedLatHash == (widenedLatHashMax - (widenedLatHashMax >> 1))) {
			return widenedLatHash | (hash & 0x5555555555555555L);
		} else {
			return (hash << 2) | (0b11);
		}
	}
	
	private static long zoomInNorth0(final long hash) {
		if ((hash & 0x2000000000000000L) != 0) {
			return 0;
		}			
		long widenedLatHash = ((hash & 0xaaaaaaaaaaaaaaaaL) << 2) | 0b10;
		long widenedLatHashMax = max(widenedLatHash);
		widenedLatHashMax |= (widenedLatHashMax >> 32);
		
		if (widenedLatHash == (widenedLatHashMax & 0xaaaaaaaaaaaaaaaaL) || widenedLatHash == (widenedLatHashMax - (widenedLatHashMax >> 1))) {
			return widenedLatHash | (hash & 0x5555555555555555L);
		} else {
			return (hash << 2) | (0b10);
		}
	}
	
	private static long zoomInSouth1(final long hash) {
		if ((hash & 0x2000000000000000L) != 0) {
			return 0;
		}			
		long widenedLatHash = (hash & 0xaaaaaaaaaaaaaaaaL) << 2;
		long widenedLatHashMax = max(widenedLatHash);
		widenedLatHashMax |= (widenedLatHashMax >> 32);
		
		if (widenedLatHash == (widenedLatHashMax & 0xaaaaaaaaaaaaaaaaL) || widenedLatHash == (widenedLatHashMax - (widenedLatHashMax >> 1))) {
			return widenedLatHash | (hash & 0x5555555555555555L);
		} else {
			return (hash << 2) | (0b01);
		}
	}
	
	private static long zoomInSouth0(final long hash) {
		if ((hash & 0x2000000000000000L) != 0) {
			return 0;
		}			
		long widenedLatHash = (hash & 0xaaaaaaaaaaaaaaaaL) << 2;
		long widenedLatHashMax = max(widenedLatHash);
		widenedLatHashMax |= (widenedLatHashMax >> 32);
		
		if (widenedLatHash == (widenedLatHashMax & 0xaaaaaaaaaaaaaaaaL) || widenedLatHash == (widenedLatHashMax - (widenedLatHashMax >> 1))) {
			return widenedLatHash | (hash & 0x5555555555555555L);
		} else {
			return hash << 2;
		}
	}
		
	
	/**
	 * Returns the geohash that makes up the northern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the northern half, returns the eastern one. 
	 * If there are four geohashes that make up the northern half, returns the far-eastern one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInNorth11(long hash) {
		return hash == 0b11 ? 0b11111 : zoomInNorth1(hash);
	}
	
	/**
	 * Returns the geohash that makes up the northern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the northern half, returns the eastern one. 
	 * If there are four geohashes that make up the northern half, returns the middle-eastern one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInNorth10(long hash) {
		return hash == 0b11 ? 0b11110 : zoomInNorth1(hash);
	}
	
	/**
	 * Returns the geohash that makes up the northern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the northern half, returns the western one. 
	 * If there are four geohashes that make up the northern half, returns the middle-western one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInNorth01(long hash) {
		return hash == 0b11 ? 0b11011 : zoomInNorth0(hash);
	}
	
	/**
	 * Returns the geohash that makes up the northern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the northern half, returns the western one. 
	 * If there are four geohashes that make up the northern half, returns the far-western one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInNorth00(long hash) {
		return hash == 0b11 ? 0b11010 : zoomInNorth0(hash);
	}
	
	/**
	 * Returns the geohash that makes up the southern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the southern half, returns the eastern one. 
	 * If there are four geohashes that make up the southern half, returns the far-eastern one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInSouth11(long hash) {
		return hash == 0b11 ? 0b11101 : zoomInSouth1(hash);
	}
	
	/**
	 * Returns the geohash that makes up the southern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the southern half, returns the eastern one. 
	 * If there are four geohashes that make up the southern half, returns the middle-eastern one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInSouth10(long hash) {
		return hash == 0b11 ? 0b11100 : zoomInSouth1(hash);
	}
	
	/**
	 * Returns the geohash that makes up the southern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the southern half, returns the western one. 
	 * If there are four geohashes that make up the southern half, returns the middle-western one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInSouth01(long hash) {
		return hash == 0b11 ? 0b11001 : zoomInSouth0(hash);
	}
	
	/**
	 * Returns the geohash that makes up the southern half of this geohash, or zero if we are at the maximum precision. 
	 * If there are two geohashes that make up the southern half, returns the western one. 
	 * If there are four geohashes that make up the southern half, returns the far-western one.
	 * <p>
	 * (Note: the only geohash that is longitudinally quadrisected instead of bisected is 0b11, the earth.)
	 */
	public static long zoomInSouth00(long hash) {
		return hash == 0b11 ? 0b11000 : zoomInSouth0(hash);
	}

	
	/**
	 * Shifts a geohash one spot to the east.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftEast(final long hash) {
		long lng = unwiden(hash);
		long max = max(lng);
        return hash & 0xaaaaaaaaaaaaaaaaL | widen((lng + 1) & max | (max - (max >> 1)));
	}
	
	/**
	 * Shifts a geohash one spot to the west.
	 * Results are undefined for invalid input geohashes.
	 */
	public static long shiftWest(final long hash) {
		long lng = unwiden(hash);
		long max = max(lng);
        return hash & 0xaaaaaaaaaaaaaaaaL | widen((lng - 1) | (max - (max >> 1)));
	}
	
	
	private static int invp(long i) {
		if (i == 0)
            return 32;
        int n = 1;
        if (i >>> 16 == 0) { n += 16; i <<= 16; }
        if (i >>> 24 == 0) { n +=  8; i <<=  8; }
        if (i >>> 28 == 0) { n +=  4; i <<=  4; }
        if (i >>> 30 == 0) { n +=  2; i <<=  2; }
        n -= i >>> 31;
        return n;
	}

	/**
	 * returns the maximum unwidened hash for same precision of this unwidened hash
	 */
	private static long max(long i) {
		i |= (i >> 1);
        i |= (i >> 2);
        i |= (i >> 4);
        i |= (i >> 8);
        i |= (i >> 16);
        return i;
	}
	
	/**
	 * puts a zero to the left of each bit in the lower 32 bits
	 */
	private static long widen(long low32) {
		  low32 |= low32 << 16; low32 &= 0x0000ffff0000ffffL;
		  low32 |= low32 << 8;  low32 &= 0x00ff00ff00ff00ffL;
		  low32 |= low32 << 4;  low32 &= 0x0f0f0f0f0f0f0f0fL;
		  low32 |= low32 << 2;  low32 &= 0x3333333333333333L;
		  low32 |= low32 << 1;  low32 &= 0x5555555555555555L;
		  return low32;
	}
	
	/**
	 * inverses the widen operation
	 */
	private static long unwiden(long wide) {
                            wide &= 0x5555555555555555L;
        wide ^= wide >> 1;  wide &= 0x3333333333333333L;
        wide ^= wide >> 2;  wide &= 0x0f0f0f0f0f0f0f0fL;
        wide ^= wide >> 4;  wide &= 0x00ff00ff00ff00ffL;
        wide ^= wide >> 8;  wide &= 0x0000ffff0000ffffL;
        wide ^= wide >> 16; wide &= 0x00000000ffffffffL;
        return wide;
	}
	
	
	
	private static final char[] base32 = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'b', 'c', 'd', 'e', 'f',
			'g', 'h', 'j', 'k', 'm', 'n', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' };

    
    private static final byte[] base32inv = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15, 16, -1, 17, 18, -1, 19, 20, -1, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};

	
	public static String toString(long geohash) {
	    char[] buf = new char[16];
	    int charPos = 16;
	    do {
	    	buf[--charPos] = base32[(int)(geohash & 31)];
	    	geohash >>>= 5;
	    } while (geohash != 0);
	    return new String(buf, charPos, (16 - charPos));
	}
	

	private static final byte BYTE_NEGATIVE_ONE = -1; //so that -1 doesn't default to int
	
	private static long valueForDigit(char digit) {
		if (digit > 'z') {
			throw new IllegalArgumentException("illegal digit: " + digit);
		}
		byte x = base32inv[digit];
		if (x == BYTE_NEGATIVE_ONE) {
			throw new IllegalArgumentException("illegal digit: " + digit);
		}
		return x;
	}
	
	public static long[] parseGeohashes(String[] strings) {
		int len = strings.length;
		long[] geohashes = new long[len];
		for (int i = 0; i < len; i++) {
			geohashes[i] = parseGeohash(strings[i]);
		}
		return geohashes;
	}
	
	public static long parseGeohash(String string) {
		char[] digits = string.toCharArray();
		int len = digits.length;
		if (len < 1 || len > 13) {
			throw new IllegalArgumentException("geohash too big: length=" + len);
		}
		
		long result = valueForDigit(digits[0]);

		if (len == 13 && result > 15L) {
			throw new IllegalArgumentException("geohash too big: " + string);
		}
		
		for (int i = 1; i < len; i++) {
			result = (result << 5) | valueForDigit(digits[i]);
		}
		if (!isValid(result)) {
			throw new IllegalArgumentException("invalid geohash: " + string);
		}
		return result;
	}

}
