#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <jpeglib.h>
#include <setjmp.h>
#include "martin.h"
#define BSIZE 28800000
#define _XOPEN_SOURCE 700

/* here are our X variables */
Display *dis;
int screen;
Window win;
GC gc;
XImage *x_image;
unsigned char *x_buffer;

/* here are our X routines declared! */
void init_x();
void close_x();
void redraw();


void disp (unsigned char *,int,int);

void usage ()
{
	printf("usage: font filename threshold [20-40 ish] star,framestcode [65A 97a]\n");
	exit (1);
}


static float notes[108]={
16.35,17.32,18.35,19.45,20.60,21.83,23.12,24.50,25.96,27.50,29.14,30.87,
32.70,34.65,36.71,38.89,41.20,43.65,46.25,49.00,51.91,55.00,58.27,61.74,
65.41,69.30,73.42,77.78,82.41,87.31,92.50,98.00,103.8,110.0,116.5,123.5,
130.8,138.6,146.8,155.6,164.8,174.6,185.0,196.0,207.7,220.0,233.1,246.9,
261.6,277.2,293.7,311.1,329.6,349.2,370.0,392.0,415.3,440.0,466.2,493.9,
523.3,554.4,587.3,622.3,659.3,698.5,740.0,784.0,830.6,880.0,932.3,987.8,
1047, 1109, 1175, 1245, 1319, 1397, 1480, 1568, 1661, 1760, 1865, 1976,
2093, 2217, 2349, 2489, 2637, 2794, 2960, 3136, 3322, 3520, 3729, 3951,
4186, 4435, 4699, 4978, 5274, 5588, 5920, 6272, 6645, 7040, 7459, 7902 };


int main(int argc,char *argv[])
{
	unsigned char *image1,*image2,*image3,*image4,*image5,*image6,*buff;
	short *wav;
	int *score,*quant;
	int i,x,y,j,k,l;
	char name[200];
	int dims[3],x_size,y_size,off;
	float mix;
  	DIR *dp;
  	struct dirent *ep;     

        image1=(unsigned char *)malloc(sizeof (char)*100*X_SIZE*Y_SIZE); // disp buffer
        image2=(unsigned char *)malloc(sizeof (char)*100*X_SIZE*Y_SIZE); // disp buffer
        image3=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image4=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // disp buffer
        image5=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // buffer
        image6=(unsigned char *)malloc(sizeof (char)*3*X_SIZE*Y_SIZE); // buffer
        score=(int *)malloc(sizeof (int)*3*X_SIZE*Y_SIZE); // buffer
        quant=(int *)malloc(sizeof (int)*3*X_SIZE*Y_SIZE); // buffer
        buff=(unsigned char *)malloc(sizeof (char)*3*10000*10000); // buffer


	init_x();
	read_JPEG_file ("./scan.jpg", image1, dims);
	printf ("X %d Y %d \n",dims[0],dims[1]);

	x_size=dims[0];y_size=dims[1];

	int r,xo,yo,r2;
	r=24;
	r2=r/2;

	// find the lines by running a rectangle over it half way down.
	int y_half,line_count,found_line;
	int lines[20],ls,le,frame,thresh;
	y_half=y_size/2;
	thresh=2000;


	while (line_count<15)
	{
	found_line=0;
	line_count=0;
	for (x=0;x<x_size;x++)
	{
		int tot;
		tot=0;
		for (y=0;y<r*2;y++)
		{
			int ip;
			int sp;
			sp=((y)*X_SIZE*3)+(x*3);
			ip=((y+r+y_half)*x_size*3)+(x*3);
			image2[sp]=image1[ip];
			image2[sp+1]=image1[ip+1];
			image2[sp+2]=image1[ip+2];
			tot+=image1[ip]; tot+=image1[ip+1]; tot+=image1[ip+2];
		}
		if (tot<thresh && !found_line){ found_line=1;ls=x;}
		if (tot>thresh && found_line){ found_line=0;le=x;lines[line_count]=(le+ls)/2;
			for (y=0;y<r*2;y++)
			{
				int sp;
				sp=((y)*X_SIZE*3)+(lines[line_count]*3);
				image2[sp]=128; image2[sp+1]=128; image2[sp+2]=255;
				sp=((y)*X_SIZE*3)+((lines[line_count]+1)*3);
				image2[sp]=128; image2[sp+1]=128; image2[sp+2]=255;
				sp=((y)*X_SIZE*3)+((lines[line_count]-1)*3);
				image2[sp]=128; image2[sp+1]=128; image2[sp+2]=255;
			}
			line_count++;
		}
	}
	disp (image2,frame,0);

	printf ("I found %d lines\n",line_count);
	thresh +=500;
	}

	//scanf("%c",name);	

	// run a square over it.
	for (y=r2;y<y_size-r2;y++)
	{
		for (x=r2;x<x_size-r2;x++)
		{
			long sum,ip;
			sum=0;
			for (yo=-r2;yo<r2;yo++)
			{
				ip=(3*x_size*(y+yo));
				for (xo=-r2;xo<r2;xo++)
				{
					sum+=(int)image1[ip+((x+xo)*3)];
					sum+=(int)image1[1+ip+((x+xo)*3)];
					sum+=(int)image1[2+ip+((x+xo)*3)];
				}
			}
			if (sum<200000){
				score[(y*x_size)+x]=1;
			}else{
				score[(y*x_size)+x]=0;
			}
		}	
	}

	//find the whys
	int found;
	found=0;
	int count;
	count=0;
	int xp[300],yp[300],my_notes[300];
	for (y=0;y<y_size;y++)
	{
		int this_line,ys,ye,yt;
		this_line=0;
		for (x=0;x<x_size;x++)
		{
			if (score[(y*x_size)+x]){ this_line=1;}
		}
		if (!found && this_line){ found=1;ys=y;}
		if (found && !this_line){ found=0;ye=y; yt=(ys+ye)/2;
			// do that line
			int x_found,xb,xe;
			x_found=0;
			this_line=0;
			for (x=0;x<x_size;x++)
			//for (x=0;x<389;x++)
			{
				if (score[(yt*x_size)+x]&& !this_line){ this_line=1;xb=x;}
				if (!score[(yt*x_size)+x] && this_line) {xe=x;
					xp[count]=(xe+xb)/2;
					yp[count]=(yt);
					this_line=0;
					printf ("Found %d %d \n",xp[count],yp[count]);
					count++;
				}
				
			}
		}
	}
	// draw lines and quantize notes.
	int max=0;
	int min=1000;
	int range,xx;
	for (x=0;x<count;x++)
	{
		if (xp[x]>max){max=xp[x];}
		if (xp[x]<min){min=xp[x];}
		for (xx=-r2/2;xx<r2/2;xx++)
		{
		for (y=-r2/2;y<r2/2;y++)
		{
			image1[((yp[x]+y)*3*x_size)+((xp[x]+xx)*3)]=128;
			image1[((yp[x]+y)*3*x_size)+((xp[x]+xx)*3)+1]=96;
			image1[((yp[x]+y)*3*x_size)+((xp[x]+xx)*3)+2]=180;
		}
		}
	}
	range=(max-min)/12;
	printf ("max %d min %d range %d\n",max,min,(max-min)/10);

	//int map[10]={11,12,14,16,17,19,21,23,24,26};
	static int map[]={0,2,4,5,7,9,11,12,14,16,17,19,21,23,24};

	for (x=0;x<count;x++)
	{
		int i,note;
		int best_min;
		best_min=10000;
		for (i=0;i<line_count;i++)
		{
			int dist;
			dist=xp[x]-lines[i];
			if (dist<0){dist=-dist;}
			if (dist<best_min){best_min=dist;my_notes[x]=map[i]+60;}
		}
		printf("x %d note %d \n",x,my_notes[x]);
	} 

	int nc,duration,rate;
	duration=60;
	rate=48000;
        wav=(short*)malloc(sizeof (short)*2*duration*rate); // wav 
	nc=0;
	float freq;
	float greq;
	double f[300];
	double g[300];
	double a[300];
	double b[300];
	float amp[300];
	for (x=0;x<300;x++){ amp[x]=0;a[x]=0;f[x]=0;b[x]=0;g[x]=0;}
	for (y=0;y<duration*rate;y++)
	{
		if (((long)y*(long)y_size)/((long)duration*(long)rate)==yp[nc])
		{
			freq=notes[my_notes[nc]];
			greq=notes[my_notes[nc]-12];
			printf ("note %d nc %d freq %f  y %d \n",my_notes[nc],nc,freq,y );
			amp[nc]=1;
			f[nc]=freq*2*M_PI/(float)rate;
			g[nc]=(greq+2)*2*M_PI/(float)rate;
			nc++;
		}
		float outl,outr;
		int i;
		outl=0;
		outr=0;
		for (i=0;i<count;i++)
		{
			outl+=sin(a[i])*amp[i];
			outr+=sin(b[i])*amp[i];
			amp[i]*=0.99999;
			a[i]+=f[i];
			b[i]+=g[i];
		}	
		wav[y*2]=20000*outl/9;
		wav[(y*2)+1]=20000*outr/9;
	}
	save_wav(wav,"test.wav", rate, 2, 2*rate*duration );
	

	int place;
	for (frame=0;frame<1800;frame++)
	{
	off=(y_size-Y_SIZE)*frame/1800;
	place=y_size*frame/1800;

	for (y=0;y<Y_SIZE;y++)
	{
		int ip,sp;
		ip=3*(y+off)*x_size;
		sp=3*y*X_SIZE;
		for (x=0;x<X_SIZE;x++)
		{
			int ipp,spp;
			spp=sp+(3*x);
			ipp=ip+(3*(x*x_size/X_SIZE));
			if (x<x_size && y+off<y_size)
			{
				if (score[((y+off-(r/2))*x_size)+(x-(r/2))]>100){ 
					image2[spp]=255;
					image2[spp+1]=0;
					image2[spp+2]=0;
				}else{
					image2[spp]=image1[ipp];
					image2[spp+1]=image1[ipp+1];
					image2[spp+2]=image1[ipp+2];
				}
			}

		}
	}
	//draw line
	int p;
	for (x=0;x<X_SIZE;x++)
	{
		p=(x*3)+(3*X_SIZE*(place-off));
		image2[p]=255; image2[p+1]=0; image2[p+2]=0;
		p=(x*3)+(3*X_SIZE*(1+place-off));
		if (p<3*X_SIZE*Y_SIZE) { image2[p]=255; image2[p+1]=0; image2[p+2]=0; }
		p=(x*3)+(3*X_SIZE*(2+place-off));
		if (p<3*X_SIZE*Y_SIZE) { image2[p]=255; image2[p+1]=0; image2[p+2]=0; }
		p=(x*3)+(3*X_SIZE*(3+place-off));
		if (p<3*X_SIZE*Y_SIZE) { image2[p]=255; image2[p+1]=0; image2[p+2]=0; }
	}


	disp (image2,frame,1);
	}

	scanf("%c",name);	



	close_x();

	exit(0);
}	

void disp (unsigned char *image2,int fram,int ab)
{
	int x,y;
	char *input;
	input=malloc(300);


       	for (y=0;y<Y_SIZE;y++)
       	{
               	int p=y*X_SIZE*3;
               	int XYP=X_SIZE*4*y;
               	for (x=0;x<X_SIZE;x++)
               	{
			int xpoint;
			int X_POINT;
			X_POINT=XYP+(4*x);
			xpoint=(x*3)+(p);

			x_buffer[X_POINT+2]=image2[xpoint];
			x_buffer[X_POINT+1]=image2[xpoint+1];
			x_buffer[X_POINT]=image2[xpoint+2];
                }
        }
	XPutImage(dis, win, gc, x_image, 0, 0, 0, 0, X_SIZE, Y_SIZE);
	sprintf(input,"./jpegs/image%04d.jpg",fram);
	if (ab){jayit(image2,X_SIZE, Y_SIZE, input);}
	free (input);
}


struct my_error_mgr {
  struct jpeg_error_mgr pub;	/* "public" fields */

  jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Always display the message. */
  /* We could postpone this until after returning, if we chose. */
  (*cinfo->err->output_message) (cinfo);

  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file (char * filename, unsigned char * dots, int * params)
{
  /* This struct contains the JPEG decompression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   */
  struct jpeg_decompress_struct cinfo;
  /* We use our private extension JPEG error handler.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct my_error_mgr jerr;
  /* More stuff */
  FILE * infile;		/* source file */
  JSAMPARRAY buffer;		/* Output row buffer */
  int row_stride;		/* physical row width in output buffer */

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  /* Step 1: allocate and initialize JPEG decompression object */

  /* We set up the normal JPEG error routines, then override error_exit. */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = my_error_exit;
  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input file, and return.
     */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
  }
  /* Now we can initialize the JPEG decompression object. */
  jpeg_create_decompress(&cinfo);

  /* Step 2: specify data source (eg, a file) */

  jpeg_stdio_src(&cinfo, infile);

  /* Step 3: read file parameters with jpeg_read_header() */

  (void) jpeg_read_header(&cinfo, TRUE);
  /* We can ignore the return value from jpeg_read_header since
   *   (a) suspension is not possible with the stdio data source, and
   *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
   * See libjpeg.txt for more info.
   */

  /* Step 5: Start decompressor */

  (void) jpeg_start_decompress(&cinfo);
  /* We can ignore the return value since suspension is not possible
   * with the stdio data source.
   */

  /* We may need to do some setup of our own at this point before reading
   * the data.  After jpeg_start_decompress() we have the correct scaled
   * output image dimensions available, as well as the output colormap
   * if we asked for color quantization.
   * In this example, we need to make an output work buffer of the right size.
   */ 
  /* JSAMPLEs per row in output buffer */
  row_stride = cinfo.output_width * cinfo.output_components;
  /* Make a one-row-high sample array that will go away when done with image */
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


  /* Step 6: while (scan lines remain to be read) */
  /*           jpeg_read_scanlines(...); */

  /* Here we use the library's state variable cinfo.output_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   */

  while (cinfo.output_scanline < cinfo.output_height) {
    /* jpeg_read_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could ask for
     * more than one scanline at a time if that's more convenient.
     */
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy (dots+(row_stride*cinfo.output_scanline),buffer[0],row_stride);
    /* Assume put_scanline_someplace wants a pointer and sample count. */
    /* put_scanline_someplace(buffer[0], row_stride); */

  }
  /* Step 7: Finish decompression */
  params[0]=cinfo.output_width;
  params[1]=cinfo.output_height;
  params[2]=cinfo.output_components;

  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  fclose(infile);

  /* And we're done! */
  return 1;
}

int jayit(unsigned char *screen,int image_width, int image_height, char *name)
{

int row_stride,ex,why,cmp,div,set;
unsigned char *image,**row_pointer,*cr,*cg,*cb;
row_pointer=(unsigned char **)malloc(1);

struct jpeg_compress_struct cinfo;
struct jpeg_error_mgr jerr;
FILE * outfile;		/* target file */
cinfo.err = jpeg_std_error(&jerr);
jpeg_create_compress(&cinfo);
if ((outfile = fopen(name, "wb")) == NULL) { 
	fprintf(stderr, "can't open file\n");
	exit(1);
}
jpeg_stdio_dest(&cinfo, outfile);
cinfo.image_width = image_width; 	/* image width and height, in pixels */
cinfo.image_height = image_height;
cinfo.input_components = 3;		/* # of color components per pixel */
cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
jpeg_set_defaults(&cinfo);
jpeg_set_quality(&cinfo,100,TRUE); /* limit to baseline-JPEG values */
jpeg_start_compress(&cinfo, TRUE);

  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */
    row_pointer[0] = & screen[cinfo.next_scanline * row_stride];
    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
jpeg_finish_compress(&cinfo);
fclose(outfile);
jpeg_destroy_compress(&cinfo);
}

void init_x()
{
/* get the colors black and white (see section for details) */
        unsigned long black,white;

        x_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //y_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        //z_buffer=(unsigned char *)malloc(sizeof(unsigned char)*4*X_SIZE*Y_SIZE);
        dis=XOpenDisplay((char *)0);
        screen=DefaultScreen(dis);
        black=BlackPixel(dis,screen),
        white=WhitePixel(dis,screen);
        win=XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,
                X_SIZE, Y_SIZE, 5, white,black);
        XSetStandardProperties(dis,win,"image","images",None,NULL,0,NULL);
        gc=XCreateGC(dis, win, 0,0);
        XSetBackground(dis,gc,black); XSetForeground(dis,gc,white);
        XClearWindow(dis, win);
        XMapRaised(dis, win);
        //XMoveWindow(dis, win,window_x,100);
        Visual *visual=DefaultVisual(dis, 0);
        x_image=XCreateImage(dis, visual, DefaultDepth(dis,DefaultScreen(dis)), ZPixmap, 0, x_buffer, X_SIZE, Y_SIZE, 32, 0);
};

void close_x() {
        XFreeGC(dis, gc);
        XDestroyWindow(dis,win);
        XCloseDisplay(dis);
        exit(1);
};

void redraw() {
        XClearWindow(dis, win);
};

