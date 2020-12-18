
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "nvutility.h"

#include "shapefil.h"
#include "version.h"


/*****************************************************************************\

    This program is public domain software that was developed by 
    the U.S. Naval Oceanographic Office.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

\*****************************************************************************/


/*

  Module Name:        read_coast

  Programmer(s):      Jan C. Depner (area.based.editor@gmail.com)

  Date Written:       July 2006

  Purpose:     build_coast reads ESRI shape files (.shp in combination with .shx) and orders them by one-degree cells
               in a bit packed, binary format.  It is a 2 pass process that is used for storing whole world coastline
               files.  Since all data is written using unsigned character buffers there are no "endian" issues to 
               deal with.  The final file is in the following format:


               Version:

                   128 characters ASCII


               Header:

                   180 X 360 groups of three, 32 bit integers stored as character buffers.  Each group consists of the
                   address of the block of segments for the associated one-degree cell, the total number of segments in
                   the cell, and the total number of vertices in the cell.  The order of the 180 X 360 cells goes from
                   west to east, south to north beginning at -90/-180.  That is, the first cell is -90/-180, the next is 
                   -90/-179, etc. until we reach -90/179 at which point we go to -89/-180 and so on.  We add 90 to all
                   latitudes and 180 to all longitudes so that we can work in positive numbers so, in essence we really
                   go from 0/0 to 89/359.


               Cell records:

                   Each cell consists of a number of segments.  Each segment in the cell consists of the following:


                       5 bits:        count bits
                       5 bits:        lon offset bits
                       5 bits:        lat offset bits
                       count bits:    count of vertices in the segment
                       18 bits:       lon bias + 2**17
                       18 bits:       lat bias + 2**17
                       26 bits:       start lon (times 100000)
                       25 bits:       start lat (times 100000)

                       count * (lon offset bits + lat offset bits):   lon and lat offsets (plus biases) from previous point


                   Some notes on the sizes - count bits is the number of bits needed to store the count of vertices in the
                   current segment.  This should never exceed 32.  If it does we've got a problem.  Lat and lon offset
                   bits is the number of bits needed to store the delta coded offsets between each lat and lon and the
                   previous lat and lon (or the start lat and lon) in the segment.  This also should never exceed 32.  The
                   lat and lon biases are values that we add to each lat and lon offset in order to make sure that we store
                   all offsets as positive values (I hate messing with sign extension ;-)  We use 18 bits to store these
                   because they should never exceed +-100000.  If they do then you have a line segment that is greater than
                   one degree and that is probably bogus.  We add 2**17 (131071) to this number before we store it in order
                   to ensure that it is always positive (stinkin' sign extension again ;-)  We store the start lat and lon
                   in 25 and 26 bits respectively because they are stored as positive integers times 100000 (range
                   0-17999999 and 0-35999999 respectively).  This gives us a resolution of about 1 meter at the equator.

                   In the olden days (when dinosaurs roamed the earth) I would have stored the start lat and lon as offsets
                   from the corner of the cell in order to save those few bits.  I would have also computed the number of
                   bits needed to store all of the bit counts.  This is what I did in 1981 with the original CIA WDBII
                   data and again in 1989 with the WVS data.  We used to have to worry about every bit back then.  Now I
                   find it much easier to understand if I keep that kind of logistical nightmare to a minimum ;-)  The
                   savings in storage are pretty minimal anyway.


  Caveats:     Requires shapelib version 1.2.10 or newer (many thanks to Frank Warmerdam for the library).


  Arguments:   Input file name(s) followed by output file name, for example:

                   build_coast gshhs_land.shp gshhs_lake.shp gshhs_isle.shp gshhs_pond.shp gshhs_all.ccl

*/


int32_t main (int32_t argc, char **argv)
{
  SHPHandle         shpHandle;
  SHPObject         *shape = NULL;
  FILE              *fp, *ofp;
  int32_t           i, j, k, m, type, numShapes, numParts, total, londeg[2], latdeg[2], diff_x[2], diff_y[2], num_vertices;
  int32_t           segCount, *segx, *segy, lastx, lasty, percent, old_percent, address, offset, xoff, yoff, num_segments;
  int32_t           range_x, range_y, count_bits, lon_offset_bits, lat_offset_bits, size, bias_x, bias_y, pos, max_bias;
  int32_t           input_file_count;
  uint8_t           start_segment = NVFalse;
  double            minBounds[4], maxBounds[4], lon[2], lat[2];
  char              fname[512], version[128], outname[512];
  uint8_t           *buffer, head_buf[12];


  printf ("\n\n%s\n\n", VERSION);


  if (argc < 3)
    {
      fprintf (stderr, "Usage: %s INPUT_FILE.shp OUTPUT_FILE\n", argv[0]);
      fprintf (stderr, "If the output file name does not have a .ccl extension it will be added.\n");
      exit (-1);
    }


  input_file_count = argc - 2;


  /*  Initialize variables  */

  total = 0;
  segx = NULL;
  segy = NULL;
  fp = NULL;
  londeg[0] = -999;
  londeg[1] = -999;
  latdeg[0] = -999;
  latdeg[1] = -999;
  lon[0] = -999.0;
  lon[1] = -999.0;
  lat[0] = -999.0;
  lat[1] = -999.0;


  /*  Make sure we don't have any old cell files hanging around since we're going to append to them  */

  for (i = 0 ; i < 180 ; i++)
    {
      for (j = 0 ; j < 360 ; j++)
        {
          sprintf (fname, "cell_%03d_%03d", j, i);
          remove (fname);
        }
    }


  for (m = 0 ; m < input_file_count ; m++)
    {
      /*  Initialize loop variables  */

      segCount = 0;
      percent = 0;
      old_percent = -1;


      /*  Open shape file  */

      shpHandle = SHPOpen (argv[m + 1], "rb");

      if (shpHandle == NULL)
        {
          perror (argv[m + 1]);
          exit (-1);
        }


      fprintf (stderr,"\n\nReading %s\n\n", argv[m + 1]);
      fflush (stderr);


      /*  Get shape file header info  */

      SHPGetInfo (shpHandle, &numShapes, &type, minBounds, maxBounds);


      /*  Read all shapes  */

      for (i = 0 ; i < numShapes ; i++)
        {
          shape = SHPReadObject (shpHandle, i);

          total += shape->nVertices;


          /*  Get all vertices  */

          if (shape->nVertices >= 2)
            {
              for (j = 0, numParts = 1 ; j < shape->nVertices ; j++)
                {
		  start_segment = NVFalse;


		  /*  Check for start of a new segment.  */

		  if (!j && shape->nParts > 0) start_segment = NVTrue;


		  /*  Check for the start of a new segment inside a larger group of points (this would be a "Ring" point).  */

		  if (numParts < shape->nParts && shape->panPartStart[numParts] == j)
		    {
		      start_segment = NVTrue;
		      numParts++;
		    }


                  /*  Bias lat and lon by 90 and 180 so that all points are positive  */

                  lon[1] = shape->padfX[j] + 180.0;
                  lat[1] = shape->padfY[j] + 90.0;


		  /*  Polygons are closed at the -180/180 boundary.  Since we're making a nice simple coastline we're going to discard points that 
		      are at exactly 180.000 or -180.000.  */

		  //if (shape->padfX[j] == -180.00000 || shape->padfX[j] == 180.00000) continue;
		  

                  /*  Damn boundary conditions!  */

                  if (lon[1] == 360.0)
		    {
		      fprintf (stderr,"%s %s %d %d %.11f %.11f\n",NVFFL,j,shape->padfX[j],shape->padfY[j]);
		      lon[1] = 180.0;
		    }


                  londeg[1] = (int32_t) lon[1];
                  latdeg[1] = (int32_t) lat[1];



                  /*  Changed cells (but not first time through)  */

                  if (latdeg[0] > -999 && (latdeg[1] != latdeg[0] || londeg[1] != londeg[0]))
                    {
                      /*  Start a new segment in a new cell  */

                      if (start_segment)
                        {
                          /*  Close last segment, don't add point, start new segment in new cell  */

                          fwrite (&segCount, sizeof (int32_t), 1, fp);

                          for (k = 0 ; k < segCount ; k++)
                            {
                              fwrite (&segx[k], sizeof (int32_t), 1, fp);
                              fwrite (&segy[k], sizeof (int32_t), 1, fp);
                            }

                          fclose (fp);

                          sprintf (fname, "cell_%03d_%03d", londeg[1], latdeg[1]);

                          if ((fp = fopen (fname, "ab")) == NULL)
                            {
                              perror (fname);
                              exit (-1);
                            }

                          segCount = 0;
                        }
                      else
                        {
                          /*  Add point to last segment, close last segment, start new segment in new cell with last point 
                              as first point in new segment  */

                          segCount++; 

                          fwrite (&segCount, sizeof (int32_t), 1, fp);

                          for (k = 0 ; k < segCount - 1 ; k++)
                            {
                              fwrite (&segx[k], sizeof (int32_t), 1, fp);
                              fwrite (&segy[k], sizeof (int32_t), 1, fp);
                            }

                          lastx = NINT (lon[0] * 100000.0);
                          lasty = NINT (lat[0] * 100000.0);
                          fwrite (&lastx, sizeof (int32_t), 1, fp);
                          fwrite (&lasty, sizeof (int32_t), 1, fp);

                          fclose (fp);


                          sprintf (fname, "cell_%03d_%03d", londeg[1], latdeg[1]);

                          if ((fp = fopen (fname, "ab")) == NULL)
                            {
                              perror (fname);
                              exit (-1);
                            }


                          segCount = 0;

                          segx = (int32_t *) realloc (segx, (segCount + 1) * sizeof (int32_t));
                          if (segx == NULL)
                            {
                              perror ("Allocating segx memory");
                              exit (-1);
                            }

                          segy = (int32_t *) realloc (segy, (segCount + 1) * sizeof (int32_t));
                          if (segy == NULL)
                            {
                              perror ("Allocating segy memory");
                              exit (-1);
                            }


                          /*  Add last point to current segment  */

                          segx[segCount] = NINT (lon[0] * 100000.0);
                          segy[segCount] = NINT (lat[0] * 100000.0);

                          segCount++;
                        }
                    }
                  else
                    {
                      /*  Start a new segment but not in a new cell  */

                      if (start_segment)
                        {
                          /*  First time through open the first file  */

                          if (latdeg[0] == -999)
                            {
                              sprintf (fname, "cell_%03d_%03d", londeg[1], latdeg[1]);

                              if ((fp = fopen (fname, "ab")) == NULL)
                                {
                                  perror (fname);
                                  exit (-1);
                                }
                            }
                          else
                            {
                              /*  Close last segment, start new segment  */

                              fwrite (&segCount, sizeof (int32_t), 1, fp);

                              for (k = 0 ; k < segCount ; k++)
                                {
                                  fwrite (&segx[k], sizeof (int32_t), 1, fp);
                                  fwrite (&segy[k], sizeof (int32_t), 1, fp);
                                }

                              segCount = 0;
                            }
                        }
                    }


                  segx = (int32_t *) realloc (segx, (segCount + 1) * sizeof (int32_t));
                  if (segx == NULL)
                    {
                      perror ("Allocating segx memory");
                      exit (-1);
                    }

                  segy = (int32_t *) realloc (segy, (segCount + 1) * sizeof (int32_t));
                  if (segy == NULL)
                    {
                      perror ("Allocating segy memory");
                      exit (-1);
                    }


                  /*  Add point to current segment  */

                  segx[segCount] = NINT (lon[1] * 100000.0);
                  segy[segCount] = NINT (lat[1] * 100000.0);

                  segCount++;


                  londeg[0] = londeg[1];
                  latdeg[0] = latdeg[1];
                  lon[0] = lon[1];
                  lat[0] = lat[1];
                }
            }


          percent = (int32_t) (((float) i / (float) numShapes) * 100.0);
          if (percent != old_percent)
            {
              fprintf (stderr, "%03d%% processed\r", percent);
              fflush (stderr);
              old_percent = percent;
            }


          SHPDestroyObject (shape);
        }


      /*  Close the last segment.  */

      fwrite (&segCount, sizeof (int32_t), 1, fp);

      for (k = 0 ; k < segCount ; k++)
        {
          fwrite (&segx[k], sizeof (int32_t), 1, fp);
          fwrite (&segy[k], sizeof (int32_t), 1, fp);
        }


      fprintf (stderr, "100%% processed\n\n");
      fprintf (stderr, "Total points processed = %d\n\n", total);
      fflush (stderr);


      SHPClose (shpHandle);


      /*  Force it to start a new segment.  */

      latdeg[0] = -888;
    }

  fclose (fp);

  if (segx != NULL) free (segx);
  if (segy != NULL) free (segy);


  percent = 0;
  old_percent = -1;
  total = 0;


  strcpy (outname, argv[argc - 1]);
  if (strcmp (&outname[strlen (outname) - 4], ".ccl")) sprintf (outname, "%s.ccl", argv[argc - 1]);


  fprintf (stderr,"\n\n%s\n\n", outname);
  fflush (stderr);


  if ((ofp = fopen (outname, "wb")) == NULL)
    {
      perror (outname);
      exit (-1);
    }


  /*  Write the header  */

  memset (version, 0, 128);
  sprintf (version, "%s\n", FILE_VERSION);
  fprintf(stderr,"%s\n",version);
  fflush (stderr);
  fwrite (version, 128, 1, ofp);


  /*  Initialize the header area  */

  for (i = 0 ; i < 180 ; i++)
    {
      for (j = 0 ; j < 360 ; j++)
        {
          offset = (i * 360 + j) * (3 * sizeof (int32_t)) + 128;

          address = 0;
          num_segments = 0;
          num_vertices = 0;

          fseek (ofp, offset, SEEK_SET);

          pos = 0;
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), address); pos += (8 * sizeof (int32_t));
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), num_segments); pos += (8 * sizeof (int32_t));
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), num_vertices);

          fwrite (head_buf, 3 * sizeof (int32_t), 1, ofp);
        }
    }


  max_bias = (int32_t) (pow (2.0, 17.0) - 1.0);


  for (i = 0 ; i < 180 ; i++)
    {
      for (j = 0 ; j < 360 ; j++)
        {
          num_segments = 0;
          num_vertices = 0;

          sprintf (fname, "cell_%03d_%03d", j, i);

          if ((fp = fopen (fname, "rb")) != NULL)
            {
              offset = (i * 360 + j) * (3 * sizeof (int32_t)) + 128;
              address = ftell (ofp);

              while (fread (&segCount, sizeof (int32_t), 1, fp))
                {
		  fprintf (stderr,"%s %s %d %d\n",NVFFL,segCount);
                  /*  Just in case we happened to write an empty segment between files ;-)  */

                  if (segCount)
                    {
                      num_vertices += segCount;
                      num_segments++;
                      total += segCount;

                      segx = (int32_t *) calloc (segCount, sizeof (int32_t));
                      if (segx == NULL)
                        {
                          perror ("Allocating segx memory");
                          exit (-1);
                        }

                      segy = (int32_t *) calloc (segCount, sizeof (int32_t));
                      if (segy == NULL)
                        {
                          perror ("Allocating segy memory");
                          exit (-1);
                        }

                      diff_x[0] = 99999999;
                      diff_x[1] = -99999999;
                      diff_y[0] = 99999999;
                      diff_y[1] = -99999999;

                      for (k = 0 ; k < segCount ; k++)
                        {
                          if (!fread (&segx[k], sizeof (int32_t), 1, fp))
			    {
			      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
			      fflush (stderr);
			      exit (-1);
			    }
                          if (!fread (&segy[k], sizeof (int32_t), 1, fp))
			    {
			      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
			      fflush (stderr);
			      exit (-1);
			    }
			  fprintf (stderr,"%s %s %d %d %d %d\n",NVFFL,k, segx[k],segy[k]);

                          if (k)
                            {
                              diff_x[0] = MIN (segx[k] - segx[k - 1], diff_x[0]);
                              diff_x[1] = MAX (segx[k] - segx[k - 1], diff_x[1]);
                              diff_y[0] = MIN (segy[k] - segy[k - 1], diff_y[0]);
                              diff_y[1] = MAX (segy[k] - segy[k - 1], diff_y[1]);
                            }
                        }


                      bias_x = -diff_x[0];
                      bias_y = -diff_y[0];

                      if (bias_x > max_bias || bias_x < -max_bias)
                        {
                          fprintf (stderr, "\n\nlon bias out of range, terminating!\n\n");
                          fprintf (stderr, "%d %d %d\n", i,j,bias_x);
                          exit (-1);
                        }


                      if (bias_y > max_bias || bias_y < -max_bias)
                        {
                          fprintf (stderr, "\n\nlat bias out of range, terminating!\n\n");
                          fprintf (stderr, "%d %d %d\n", i,j,bias_y);
                          exit (-1);
                        }


                      range_x = diff_x[1] - diff_x[0];
                      range_y = diff_y[1] - diff_y[0];


                      if (!range_x) range_x = 1;
                      if (!range_y) range_y = 1;


                      count_bits = NINT (log10 ((double) segCount) / log10 (2.0) + 1.0);
                      lon_offset_bits = NINT (log10 ((double) range_x) / log10 (2.0) + 1.0);
                      lat_offset_bits = NINT (log10 ((double) range_y) / log10 (2.0) + 1.0);


                      size = 5 + 5 + 5 + count_bits + lon_offset_bits + lat_offset_bits + 18 + 18 + 26 + 25 + 
                        (segCount - 1) * (lon_offset_bits + lat_offset_bits);

                      size = size / 8 + 1;

                      buffer = (uint8_t *) calloc (1, size);

                      if (buffer == NULL)
                        {
                          perror ("Allocating buffer");
                          exit (-1);
                        }

                      pos = 0;
                      bit_pack (buffer, pos, 5, count_bits); pos += 5;
                      bit_pack (buffer, pos, 5, lon_offset_bits); pos += 5;
                      bit_pack (buffer, pos, 5, lat_offset_bits); pos +=5;
                      bit_pack (buffer, pos, count_bits, segCount); pos += count_bits;
                      bit_pack (buffer, pos, 18, bias_x + max_bias); pos += 18;
                      bit_pack (buffer, pos, 18, bias_y + max_bias); pos += 18;
                      bit_pack (buffer, pos, 26, segx[0]); pos += 26;
                      bit_pack (buffer, pos, 25, segy[0]); pos += 25;


                      for (k = 1 ; k < segCount ; k++)
                        {
                          xoff = (segx[k] - segx[k - 1]) + bias_x;
                          yoff = (segy[k] - segy[k - 1]) + bias_y;

                          bit_pack (buffer, pos, lon_offset_bits, xoff); pos += lon_offset_bits;
                          bit_pack (buffer, pos, lat_offset_bits, yoff); pos += lat_offset_bits;
                        }

                      fwrite (buffer, size, 1, ofp);

                      free (buffer);
                      free (segx);
                      free (segy);
                    }
                }

              fclose (fp);
              remove (fname);


              /*  Write the address, number of segments, and number of vertices in the header  */

              fseek (ofp, offset, SEEK_SET);

              k = 8 * sizeof (int32_t);

              pos = 0;
              bit_pack (head_buf, pos, k, address); pos += k;
              bit_pack (head_buf, pos, k, num_segments); pos += k;
              bit_pack (head_buf, pos, k, num_vertices);

              fwrite (head_buf, 3 * sizeof (int32_t), 1, ofp);

              fseek (ofp, 0, SEEK_END);
            }
        }

      percent = (int32_t) (((float) i / 181.0) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "%03d%% packed\r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  fclose (ofp);

  fprintf (stderr, "100%% packed\n\n");
  fprintf (stderr, "Total points packed = %d\n\n", total);
  fflush (stderr);

  return (0);
}
