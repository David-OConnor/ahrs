use core::f32::consts::TAU;

use num_traits::Float;

// todo: Update the table; it's from 2012

// #include <stdio.h>
// #include <math.h>

// /*
// 	Adjust electronic compass readings from magnetic north to true north using a subset of data generated by the NOAA WMM70.exe software.
// 	Approximate magnetic declination for latitudes -60..60 and all longitudes for July 15, 2012 (the midpoint of the the period covered by
// 	the model.) Points at 10 degree intervals of a (latitude, longitude) grid are plotted for declination	with the WMM2010 software.
// 	Declinations are stored in the source program in a literal string. Bilinear interpolation is then used to approximate the declinations
// 	within these 10 degree grids. This method allows microcontrollers with minimal memory resources to correct readings from a compass sensor
// 	for true north. Test readings should be within a few degrees of the WMM2010 model for most inhabited locations.

// 	If you just need to use the existing lookup table here for use in your program, you will need
// 	only the	declination() function to retrieve values from the lookup table dec_tbl[]
// 	Call declination with your latitude and longitude. approximation for declination is returned.

// 	To create your own lookup table:

// 	STEP 1: 	download and install geomag70.exe from 	www.ngdc.noaa.gov/geomag/WMM/soft.shtml

// 	STEP 2: 	read the readme.txt file
// 				eg. run the program at the command prompt for help: geomag70 h
// 				then modify fprintf() format string in the build_wmm_input_file() for your required values
// 				run this program and choose '3. Build Input File' to create geomag70_input.txt
// 				the file contents will look something like this:

// 				2012.5 D M100 -60 -180
// 				2012.5 D M100 -60 -170
// 				:
// 				2012.5 D M100 60 170
// 				2012.5 D M100 60 180

// 	STEP 3: 	To create the raw output file run: geomage70.exe  WMM2010.COF f geomag70_input.txt geomag70_output.txt
// 				the file contents will look something like this:

// 				Date Coord-System Altitude Latitude Longitude D_deg D_min I_deg I_min H_nT X_nT Y_nT Z_nT F_nT dD_min dI_min dH_nT dX_nT dY_nT dZ_nT dF_nT
// 				2012.5 D M100 -60 -180    46d 37m   -77d 45m   13243.3   9095.8   9625.6 -60953.3  62375.4     6.6       1.2         10.5    -11.2     25.0     55.3    -51.8
// 				2012.5 D M100 -60 -170    45d 31m   -75d 45m   14906.9  10446.2  10634.5 -58725.1  60587.5     5.8       1.4          8.1    -12.3     23.3     66.1    -62.1
// 				:
// 				2012.5 D M100 60 170    -1d 28m    70d 21m   18169.1  18163.1   -467.1  50891.8  54037.8    -2.2       1.9         -6.5     -6.8    -11.6     68.4     62.3
// 				2012.5 D M100 60 180     3d 58m    70d 23m   18041.0  17997.9   1247.3  50595.9  53716.1    -6.4       1.8         -5.0     -2.7    -33.6     69.0     63.3

// 	STEP 4:	To build the lookup table run this program and choose '4. Build Lookup Table'
// 				this will extract the declination values from the 6th and 7th columns in the geomag70_output.txt file for each lat,lon
// 				to create the file wmm_ss.txt containing a comma-separated list of declinations as minutes
// 				the file contents will look something like this:

// 				46,45,44,42,41,40,37,-68,0,0,0,...

// 				Break this into multiple lines in your program with the '\' line continuation directive:

// 				46,45,...32,-4, \    (one space after comma, then '\', then hit return)
// 				37,22-3,...,25, \
// 				etc.

// 	STEP 5:	using a text editor, copy the data into the declination() function dec_tbl[] lookup table

// 	STEP 6: 	test the declination() function by choosing option '2. Lookup Declinations' using various latitudes and longitudes
// 				compare the approximated declinations with the geomag_70.exe calculated declinations

// 	scottfromscott@bellsouth.net
// */
// void build_wmm_input_file()
// {
// 	/*
// 		Build input file for processing by the geomag70.exe program available from the NOAA web site above
// 		The input file has the following format
// 		where columns are date, coord_system, meters_above_sea_level, lat, lon
// 		2012.5 is midpoint of WMM2010 model's range
// 		100 meters above sea level is is the average elevation of populated areas

// 		2012.5 D M100 -60 -180
// 		2012.5 D M100 -60 -170
// 		:
// 		2012.5 D M100 60 170
// 		2012.5 D M100 60 180
// 	*/
// 	int lat_index, lon_index;

// 	FILE *fp;

// 	fp = fopen("geomag70_input.txt", "w");
// 	for (lat_index=-60;lat_index<=60;lat_index += 10)
// 	{
// 		for(lon_index=-180;lon_index<=180;lon_index += 10)
// 		{
// 			fprintf(fp,"2012.5 D M100 %hd %hd\n",lat_index,lon_index);
// 		}
// 	}
// 	fclose(fp);
// }

// void build_lookup_table()
// {
// 	/*
// 		for each lat,lon, read declination value from geomag70_output.txt file
// 		write declination value rounded to nearest degree as char type: -128..127 to the wmm_ss.txt file
// 	*/
// 	FILE *inf, *outf;

// 	short d, m, dec, entry_counter;
// 	int i, j, k, lat_index, lon_index;
// 	char line[160], dec_deg[4], dec_min[4];

// 	inf = fopen("geomag70_output.txt", "r");
// 	outf = fopen("wmm_ss.txt", "w");

// 	fgets ( line, sizeof(line), inf );

// 	// -60..0..60 step 10 => 2(6) + 1 = 13 dimensions for lat
// 	for (lat_index=0;lat_index<13;lat_index++)
// 	{// -180..0..180 step 10 => 2(18) + 1 = 37 dimensions for lon
// 		for(lon_index=0;lon_index<37;lon_index++)
// 		{
// 			// get declination from file

// 			fgets ( line, sizeof(line), inf );

// 			// dig out the degrees and minutes

// 			i = 0;
// 			for(j=1;j<=5;j++)
// 			{
// 				while(line[i]!=' ') i++;
// 				i++;
// 			}
// 			i+=2;
// 			j = 0;
// 			while (line[i]!='d') dec_deg[j++] = line[i++];
// 			dec_deg[j]='\0';
// 			i+=2;
// 			j=0;
// 			while (line[i]!='m') dec_min[j++] = line[i++];
// 			dec_min[j]='\0';

// 			d = atoi(dec_deg); // d is signed
// 			m = atoi(dec_min);
// 			m=(d<0)? -m : m; // m is not signed, so adjust it for negative d's
// 			dec = round(d+(m/60)); // round to nearest degree

// 			//output to clean file...
// 			fprintf(outf,"%i,",dec);
// 		}
// 	}
// 	fclose(inf);
// 	fclose(outf);
// }

// -60..0..60 step 10 => 2(6) + 1 = 13 dimensions for lat; -180..0..180 step 10 => 2(18) + 1 = 37 dimensions for lon
#[rustfmt::skip]
const LUT: [i8; 13 * 37] = [
    46,45,44,42,41,40,38,36,33,28,23,16,10,4,-1,-5,-9,-14,-19,-26,-33,-40,-48,-55,-61, 
        -66,-71,-74,-75,-72,-61,-25,22,40,45,47,46,30,30,30,30,29,29,29,29,27,24,18,11,3, 
        -3,-9,-12,-15,-17,-21,-26,-32,-39,-45,-51,-55,-57,-56,-53,-44,-31,-14,0,13,21,26, 
        29,30,21,22,22,22,22,22,22,22,21,18,13,5,-3,-11,-17,-20,-21,-22,-23,-25,-29,-35, 
        -40,-44,-45,-44,-40,-32,-22,-12,-3,3,9,14,18,20,21,16,17,17,17,17,17,16,16,16,13, 
        8,0,-9,-16,-21,-24,-25,-25,-23,-20,-21,-24,-28,-31,-31,-29,-24,-17,-9,-3,0,4,7, 
        10,13,15,16,12,13,13,13,13,13,12,12,11,9,3,-4,-12,-19,-23,-24,-24,-22,-17,-12,-9, 
        -10,-13,-17,-18,-16,-13,-8,-3,0,1,3,6,8,10,12,12,10,10,10,10,10,10,10,9,9,6,0,-6, 
        -14,-20,-22,-22,-19,-15,-10,-6,-2,-2,-4,-7,-8,-8,-7,-4,0,1,1,2,4,6,8,10,10,9,9,9, 
        9,9,9,8,8,7,4,-1,-8,-15,-19,-20,-18,-14,-9,-5,-2,0,1,0,-2,-3,-4,-3,-2,0,0,0,1,3,5, 
        7,8,9,8,8,8,9,9,9,8,8,6,2,-3,-9,-15,-18,-17,-14,-10,-6,-2,0,1,2,2,0,-1,-1,-2,-1,0, 
        0,0,0,1,3,5,7,8,8,9,9,10,10,10,10,8,5,0,-5,-11,-15,-16,-15,-12,-8,-4,-1,0,2,3,2,1,0, 
        0,0,0,0,-1,-2,-2,-1,0,3,6,8,6,9,10,11,12,12,11,9,5,0,-7,-12,-15,-15,-13,-10,-7,-3, 
        0,1,2,3,3,3,2,1,0,0,-1,-3,-4,-5,-5,-2,0,3,6,5,8,11,13,15,15,14,11,5,-1,-9,-14,-17, 
        -16,-14,-11,-7,-3,0,1,3,4,5,5,5,4,3,1,-1,-4,-7,-8,-8,-6,-2,1,5,4,8,12,15,17,18,16, 
        12,5,-3,-12,-18,-20,-19,-16,-13,-8,-4,-1,1,4,6,8,9,9,9,7,3,-1,-6,-10,-12,-11,-9,-5, 
        0,4,3,9,14,17,20,21,19,14,4,-8,-19,-25,-26,-25,-21,-17,-12,-7,-2,1,5,9,13,15,16,16, 
        13,7,0,-7,-12,-15,-14,-11,-6,-1,3
];

/// lat and lon input here are in degrees. The result is in radians.
pub fn estimate_declination(lat: f32, lon: f32) -> f32 {
    /*
        lookup declination values from 10 x 10 degree grid and return approximate declination for (lat,lon)
        data; 482 declination values (rounded to nearest degree) stored in 482 bytes
    */

    /*
        DIAGRAM OF 10X10 GRID SQUARE:

    (+60)						(lon,latmin+10,decmax=?)
    l	(lonmin,latmin+10,decNW)|	|		|(lonmin+10,latmin+10,decNE)
    a								 --o--x-----o--
    t									|	|		|
    i									|	|		|
    t									+--x(lon,lat,dec=?)
    u									|	|		|
    d								 --o--x-----o--
    e		(lonmin,latmin,decSW)|	|		|(lonmin+10,latmin,decSE)
    (-60)						(lon,latmin,decmin=?)

                             (-180)<- longitude ->(+180)

                        o = decs from table, x = calculated decs
    */

    /* set base point (latmin, lonmin) of grid */

    /* no limits test on lon */
    let lonmin = if lon == 180. {
        170.
    } else {
        (lon / 10.).floor() * 10.
    };

    /* supported lat's -60..60, so range check... */
    let latmin = if lat >= 60. {
        60.
    } else if lat < -60. {
        -60.
    } else {
        (lat / 10.).floor() * 10.
    };

    /* array index = (degrees+[60|180])/10 */
    let latmin_index = ((60. + latmin) as usize) / 10;
    let lonmin_index = ((180. + lonmin) as usize) / 10;

    // todo: QC row/col order
    let dec_sw = LUT[latmin_index * 37 + lonmin_index] as f32;
    let dec_se = LUT[latmin_index * 37 + lonmin_index + 1] as f32;
    let dec_ne = LUT[(latmin_index + 1) * 37 + lonmin_index + 1] as f32;
    let dec_nw = LUT[(latmin_index + 1) * 37 + lonmin_index] as f32;

    /* approximate declination within the grid using bilinear interpolation */

    let decmin = (lon - lonmin) / 10. * (dec_se - dec_sw) + dec_sw;
    let decmax = (lon - lonmin) / 10. * (dec_ne - dec_nw) + dec_nw;

    let result_deg = (lat - latmin) / 10. * (decmax - decmin) + decmin;
    result_deg * TAU / 360.
}
