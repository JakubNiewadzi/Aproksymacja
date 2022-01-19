#ifndef MAKESPL_H
#define MAKESPL_H

#include "points.h"
#include "splines.h"


void  make_spl ( points_t *pts, spline_t *spl);
void our_make(points_t* pts, char* out, double fromX, double toX, int n);

#endif
