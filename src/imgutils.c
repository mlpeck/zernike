/* jpeg, tiff file input and basic image processing routines.
  jpeg and tiff input functions adapted from package biOps,
  Copyright 2007 Walter Alini, Matías Bordese.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* get lineal pos from x,y,c(hannel),w(idth),h(eight) */
#define IMGPOS(x,y,c,w,h) ((h)*(x) + (y) + (w)*(h)*(c))
#define MATPOS(x,y,h) ((h)*(x) + (y))

/* min and max macros */

#define MIN(a, b) (a) >= (b) ? (b) : (a)
#define MAX(a, b) (a) >= (b) ? (a) : (b)

void resize_image(double *im_in, int *width, int *height, 
  double *im_out, int *w_out, int *h_out, int *ret) {
  
  double si_w = (double) ((*width)-1) / (double) ((*w_out)-1);
  double si_h = (double) ((*height)-1) / (double) ((*h_out)-1);
  
  int i,j, ii0, ii1, jj0, jj1;
  double p, q, i_old, j_old;
  
  for (i=0; i<(*w_out); i++) {
	i_old = si_w*i;
	ii0 = floor(i_old);
	ii1 = ii0+1;
	p = i_old - (double) ii0;
	for (j=0; j<(*h_out); j++) {
	  j_old = si_h*j;
	  jj0 = floor(j_old);
	  jj1 = jj0+1;
	  q = j_old - (double) jj0;
	  im_out[MATPOS(i,j,(*h_out))] = 
		(1.0-p)*(1.0-q)*im_in[MATPOS(ii0,jj0,(*height))] +
		p*(1.0-q)*im_in[MATPOS(MIN(ii1,(*width-1)),jj0,(*height))] +
		(1.0-p)*q*im_in[MATPOS(ii0,MIN(jj1,(*height-1)),(*height))] +
		p*q*im_in[MATPOS(MIN(ii1,(*width-1)),MIN(jj1,(*height-1)),(*height))];
	}
  }
  *ret=1;
}

void comb_channels(int *img_in, double *img_out, 
		int *width, int *height, int *depth, double *channels) {
	
	int i, j, k;
	double imgval;
	
	for (i=0; i < (*width); i++) {
	  for (j=0; j < (*height); j++) {
		imgval = 0.0;
		for (k=0; k < (*depth); k++)
		  imgval += channels[k] * img_in[IMGPOS(i, j, k, (*width), (*height))];
		img_out[MATPOS(i,j,(*height))] = imgval;
	  }
	}
}

/* Copyright 2007 Walter Alini, Matías Bordese */

#include <jpeglib.h>

/*
	Function: read_jpg_img_info
	Read the image info: width, height and depth.

	Parameters:
		filename - The file path of the image.

	Returns:
		width - The width of the image.
		height - Its height.
		depth - Its color depth.
		ret - -1 if there is an error, 1 if success.
*/
void read_jpg_img_info(char **filename, int *width, int *height, int *depth, int *ret){
	char *fname = *filename;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;
	
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	
	if ((infile = fopen(fname, "rb")) == NULL) {
		*ret = -1;    /* couldn't open file */
		return;
	}

	jpeg_stdio_src(&cinfo, infile);
	jpeg_read_header(&cinfo, TRUE);

	*width = cinfo.image_width;
	*height = cinfo.image_height;
	*depth = cinfo.num_components;

	jpeg_destroy_decompress(&cinfo);
	fclose(infile);
	*ret = 1;
}

/*
	Image data format:
	[[R] [G] [B]] where [X] = [[col_1] ... [col_width]] where [col_i] = [row_1 ... row_height]
*/

/*
	Function: read_jpg_img
	Read the image data.

	Parameters:
		filename - The file path of the image.

	Returns:
		image - The image data.
		ret - -1 if there is an error, 1 if success.
*/
void read_jpg_img(char **filename, int *width, int *height, int *depth, 
				  double *channels, double *image, int *ret){
	char *fname = *filename;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;		/* source file */
	JSAMPARRAY buffer;	/* Output row buffer */
	int row_stride;		/* physical row width in output buffer */
	int plane_size;
	int i, j;
	int line;
	unsigned char *p;
	int *img_in = malloc((*width)*(*height)*(*depth)*sizeof(int));
	
	if ((infile = fopen(fname, "rb")) == NULL) {
		*ret = -1;   /* couldn't open file */
		return;
	}

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, infile);

	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);
	row_stride = cinfo.output_width * cinfo.output_components;
	plane_size = cinfo.output_width * cinfo.output_height;
	buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
	
	while (cinfo.output_scanline < cinfo.output_height) {
		line = cinfo.output_scanline;  /* preserve current scanline */
		jpeg_read_scanlines(&cinfo, buffer, 1);
		p = buffer[0];
		for (i = 0; i < cinfo.output_width; i++) {
			for (j = 0; j < cinfo.output_components; j++) {
				img_in[line + cinfo.output_height * i + j * plane_size] = (unsigned char) *p++;
			}
		}
	}

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(infile);
	if ((*depth) == 1) {
	  for (i=0; i<(*width) * (*height); i++) image[i] = (double) img_in[i];
	}
	  else comb_channels(img_in, image, width, height, depth, channels);
	*ret = 1;
    free(img_in);
	return;
}

#ifdef TIFF_OK	/* Kludge until I figure out how to link zlib correctly */

#include <tiffio.h>

/*
	Function: read_tiff_img_info
	Read the image info: width, height and depth.

	Parameters:
		filename - The file path of the image.

	Returns:
		width - The width of the image.
		height - Its height.
		depth - Its color depth.
		ret - -1 if there is an error, 1 if success.
*/
void read_tiff_img_info(char **filename, int *width, int *height, int *depth, int *ret){
	char *fname = *filename;
	TIFF *image;
	
	if ((image = TIFFOpen(fname, "r")) == NULL) {
		*ret = -1;    /* couldn't open file */
		return;
	}

	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, width);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, height);
	*depth = 3;

	TIFFClose(image);
	*ret = 1;
}

/*
	Image data format:
	[[R] [G] [B]] where [X] = [[col_1] ... [col_width]] where [col_i] = [row_1 ... row_height]
*/

/*
	Function: read_tiff_img
	Read the image data.

	Parameters:
		filename - The file path of the image.

	Returns:
		image - The image data.
		ret - -1 if there is an error, 1 if success.
*/
void read_tiff_img (char **filename, int *width, int *height, int *depth,
					double *channels, double *image, int *ret){
	char *fname = *filename;
	TIFF *tiff_image;
	uint32 twidth, theight, *raster;
	unsigned long imagesize, i, j, line, plane_size;
	int *img_in = malloc((*width)*(*height)*(*depth)*sizeof(int));

	if ((tiff_image = TIFFOpen(fname, "r")) == NULL) {
		*ret = -1;    /* couldn't open file */
		return;
	}

	/* Find the width and height of the image */
	TIFFGetField(tiff_image, TIFFTAG_IMAGEWIDTH, &twidth);
	TIFFGetField(tiff_image, TIFFTAG_IMAGELENGTH, &theight);
	imagesize = theight * twidth + 1;

	if ((raster = (uint32 *) malloc(sizeof(uint32) * imagesize)) == NULL){
		fprintf(stderr, "Could not allocate enough memory\n");
		*ret = -1;
		return;
	}

	/* Read the image into the memory buffer */
	if (TIFFReadRGBAImage(tiff_image, twidth, theight, raster, 0) == 0){
		fprintf(stderr, "Could not read image\n");
		*ret = -1;
		return;
	}

	line = 0;
	plane_size = twidth * theight;
	for (i = theight - 1; i != -1; i--){
		for (j = 0; j < twidth; j++){
			img_in[line + theight * j] = (unsigned char) TIFFGetR(raster[i * twidth + j]);
			img_in[line + theight * j + plane_size] = (unsigned char) TIFFGetG(raster[i * twidth + j]);
			img_in[line + theight * j + plane_size * 2] = (unsigned char) TIFFGetB(raster[i * twidth + j]);
		}
		line++;
	}

	free(raster);
	TIFFClose(tiff_image);
	if ((*depth) == 1) {
	  for (i=0; i<(*width) * (*height); i++) image[i] = (double) img_in[i];
	}
	  else comb_channels(img_in, image, width, height, depth, channels);

	*ret = 1;
	return;
}

#endif
