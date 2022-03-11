//
// Created by ac126718 on 11.03.2022.
//

#include "matplotHelper.h"

namespace Ikarus::plot
{
  void draw_xy(const Eigen::VectorXd&x,const Eigen::VectorXd& y)
  {
    assert(x.size()==y.size() && "The passed x and y vectors have to have the same size!");
    std::vector<double> xstd(x.begin(),x.end());
    std::vector<double> ystd(y.begin(),y.end());
    matplot::plot(x, y, "-o");
    matplot::hold(matplot::on);
  }


}