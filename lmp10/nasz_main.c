#include "points.h"
#include "makespl.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

char *usage =
  "Usage: %s  -p points-file [ -g gnuplot-file [-f from_x -t to_x -n n_points ] ]\n"
  "            if points-file is given then\n"
  "               reads discrete 2D points from points-file\n"
  "               - number of points should be >= 4\n"
  "            else (points-file not given)\n"
  "            if gnuplot-file is given then\n"
  "               makes table of n_points within <from_x,to_x> range\n"
  "               - from_x defaults to x-coordinate of the first point in points-file,\n"
  "               - to_x defaults to x-coordinate of the last point\n"
  "               - n_points defaults to 100\n"
  "               - n_points must be > 1\n"
  "            endif\n";


int
main (int argc, char **argv)
{
  int opt;
  char *inp = NULL;
  char *out = NULL;
  char *gpt = NULL;
  double fromX = 0;
  double toX = 0;
  int n = 0;
	char *progname= argv[0];

  points_t pts;

  pts.n = 0;
  /* process options, save user choices */
  while ((opt = getopt (argc, argv, "p:g:f:t:n:")) != -1) {
    switch (opt) {
    case 'p':
      inp = optarg;
      break;
    case 'g':
      gpt = optarg;
      break;
    case 'f':
      fromX = atof (optarg);
      break;
    case 't':
      toX = atof (optarg);
      break;
    case 'n':
      n = atoi (optarg);
      break;
    default:                   /* '?' */
      fprintf (stderr, usage, progname);
      exit (EXIT_FAILURE);
    }
  }
	if( optind < argc ) {
		fprintf( stderr, "\nBad parameters!\n" );
		for( ; optind < argc; optind++ )
			fprintf( stderr, "\t\"%s\"\n", argv[optind] );
		fprintf( stderr, "\n" );
		fprintf( stderr, usage, progname );
		exit( EXIT_FAILURE );
	}
    if (inp != NULL) {
    FILE *ouf = NULL; 

    FILE *inf = fopen (inp, "r");
    if (inf == NULL) {
      fprintf (stderr, "%s: can not read points file: %s\n\n", argv[0], inp);
      exit (EXIT_FAILURE);
    }

    if (read_pts_failed (inf, &pts)) {
      fprintf (stderr, "%s: bad contents of points file: %s\n\n", argv[0],
               inp);
      exit (EXIT_FAILURE);
    }
    else
      fclose (inf);
}else{
  exit(EXIT_FAILURE);
}
printf("%lf, %lf, %d %s", fromX, toX, n, gpt);


  our_make(&pts, gpt, fromX, toX, n);
  return 0;
  
}
