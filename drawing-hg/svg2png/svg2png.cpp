#include <QtGui>
#include <QSvgRenderer>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

int main(int argc, char* argv[])
{
  if (argc < 5) {
    printf("Usage %s <width> <height> <infile> <outfile>\n", argv[0]);
    return 1;
  }
  QApplication app(argc, argv, false);
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);

  // protection against attempts to pass too large values
  if (width > 8192 && height > 8192) {
    fprintf(stderr, "Image too large.\n");
    return 1;
  }

  QSvgRenderer renderer(QString::fromLocal8Bit(argv[3]));
  if (height==0){
    QSize size=renderer.defaultSize ();
    height = (int)(((float)size.height())/((float)size.width())*(float)width);
  }

  if (width <= 0|| height <= 0) {
    fprintf(stderr, "Please specify a valid size.\n");
    return 1;
  }


  QImage img(width, height, QImage::Format_ARGB32);
  img.fill(0);

  if (!renderer.isValid()){
    std::cerr<<"invalid document detected!"<<std::endl;
    return 0;
  }

  QPainter p(&img);
  p.setRenderHint(QPainter::Antialiasing, true);
  p.setRenderHint(QPainter::TextAntialiasing, true);
  p.setRenderHint(QPainter::SmoothPixmapTransform, true);

  renderer.render(&p);

  img.save(QString::fromLocal8Bit(argv[4]), "PNG");
  return 0;
}
