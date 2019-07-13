
#ifndef DIST_HEATMAP_H
#define DIST_HEATMAP_H


#include "sam.h"
#include "paris.h"
#include "paris_plot.h"
#include <string_split.h>

using DistMatrix = Matrix<double>;

void read_matrix(const string &file_name, DistMatrix &dist_matrix);

void compress_matrix_to_heatmap(const DistMatrix &dist_matrix,
                                pan::Heatmap_Matrix &heatmap_matrix,
                                uLONG target_size, 
                                double dist_1,
                                double dist_2,
                                double dist_3,
                                string color_1,
                                string color_2,
                                string color_3,
                                string color_4);

#endif