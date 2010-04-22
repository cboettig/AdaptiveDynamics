#include "png_plots.hh"

void boxplot(int x, int y, pngwriter &png)
{
	png.plot(x+1,y, 0., 0., 1.);
	png.plot(x-1,y, 0., 0., 1.);
	png.plot(x,y+1, 0., 0., 1.);
	png.plot(x,y-1, 0., 0., 1.);

	png.plot(x-1,y-1, 0., 0., 1.);
	png.plot(x+1,y+1, 0., 0., 1.);
	png.plot(x-1,y+1, 0., 0., 1.);
	png.plot(x+1,y-1, 0., 0., 1.);
/*
	png.plot(x+2,y, 0., 0., 1.);
	png.plot(x-2,y, 0., 0., 1.);
	png.plot(x,y+2, 0., 0., 1.);
	png.plot(x,y-2, 0., 0., 1.);
*/
}


void plotpng(vector<pop> &poplist, int time, pngwriter &png)
{
	vector<pop>::iterator p;
	int x=time, y;
	double ymax = .75;
	int shift = Y_PLOTSIZE/2;

	for(p=poplist.begin(); p != poplist.end(); p++){
		x = time;
		y = shift + (int) round((p->trait/ymax)*Y_PLOTSIZE/2);
		png.plot(x, y, 0.0, 0.0, 1.0);
		boxplot(x, y, png);
	}
	/* crude gridlines */
	png.plot(x, shift, .5, .5, .5);
	double secondline = .25;
	int grid = (int) ((secondline/ymax)*Y_PLOTSIZE/2);
	png.plot(x, shift+grid, .5, .5, .5);
	png.plot(x, shift-grid, .5, .5, .5);
}


