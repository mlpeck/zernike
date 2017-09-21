# include <libraw/libraw.h>
# include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix readraw(CharacterVector fname, NumericVector channels) {
  LibRaw rawimage;
  int ret=0;
  int i, j;
  
  std::string fn = as<std::string>(fname);
  ret=rawimage.open_file(fn.c_str());
  if (ret != 0) stop("LibRaw returned error");
  ret = rawimage.unpack();
  ret = rawimage.subtract_black();
  int h = rawimage.imgdata.sizes.height;
  int w = rawimage.imgdata.sizes.width;
  
  // set gamma slopes to 1 and color space to raw
  
  rawimage.imgdata.params.gamm[0] = 1.0;
  rawimage.imgdata.params.gamm[1] = 1.0;
  rawimage.imgdata.params.output_color = 0;
  ret = rawimage.dcraw_process();
  if (ret != 0) stop("LibRaw returned error");
  NumericMatrix img(w, h);
  channels = channels/sum(channels);
  for (int k=0; k<h*w; k++) {
    i = k % w;
    j = k / w;
    img(i, j) = 0;
    for (int l=0; l<3; l++) {
      img(i, j) += channels(l) * (double) rawimage.imgdata.image[k][l];
    }
  }
  return img;  
}
