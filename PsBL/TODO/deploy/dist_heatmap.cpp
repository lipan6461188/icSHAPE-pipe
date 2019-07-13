
#include "dist_heatmap.h"

using namespace std;

void read_matrix(const string &file_name, DistMatrix &dist_matrix)
{
    dist_matrix.clear();
    ifstream IN(file_name, ifstream::in);
    string cur_line;
    while(getline(IN, cur_line))
    {
        trim(cur_line, '\t');
        auto data = split(cur_line, '\t');
        DoubleArray dist_line;
        for(string item: data)
        {
            if(item == "NULL")
                dist_line.push_back(-1);
            else{
                dist_line.push_back( stod(item) );
            }
        }
        dist_matrix.push_back(dist_line);
    }
    IN.close();

    // check
    for(const auto &row: dist_matrix)
        if(row.size() != dist_matrix.size())
            throw runtime_error(file_name+" matrix dimesion is invalid");
}


double matrix_block_mean(const DistMatrix &matrix,
                            DistMatrix::size_type x,
                            DistMatrix::size_type y,
                            DistMatrix::size_type x_len,
                            DistMatrix::size_type y_len,
                            double exclude_v = -1)
{
    double total_v = 0;
    uLONG valid_count = 0;

    for(decltype(matrix.size()) idx=x; idx<x+x_len; idx++)
        for(decltype(matrix.size()) idy=y; idy<y+y_len; idy++)
        {
            auto cur_v = matrix.at(idx).at(idy);
            if(cur_v != exclude_v)
            {
                total_v += cur_v;
                ++valid_count;
            }
        }

    if( valid_count > (x_len*y_len)/2 )
        return total_v / valid_count;
    else
        return -1;
}

void compress_matrix_to_heatmap(const DistMatrix &dist_matrix,
                                pan::Heatmap_Matrix &heatmap_matrix,
                                uLONG target_size, 
                                double dist_1,
                                double dist_2,
                                double dist_3,
                                string color_1,
                                string color_2,
                                string color_3,
                                string color_4)
{
    using size_type = DistMatrix::size_type;

    if(target_size*2 >= dist_matrix.size())
    {
        throw runtime_error("Invalid target_size to compress matrix");
    }
    
    init_matrix(heatmap_matrix, target_size);

    //QColor c_1( QString::fromStdString(color_1) ), c_2( QString::fromStdString(color_2) ), c_3( QString::fromStdString(color_3) );
    //QColor black(Qt::);

    double step = 1.0 * dist_matrix.size() / target_size;
    for(size_type idx=0; idx<target_size; idx++)
        for(size_type idy=0; idy<target_size; idy++)
        {
            size_type x_base = idx * step;
            size_type y_base = idy * step;
            size_type x_len = min( decltype(dist_matrix.size())((idx+1) * step), dist_matrix.size() ) - x_base;
            size_type y_len = min( decltype(dist_matrix.size())((idy+1) * step), dist_matrix.size() ) - y_base;
            double value = matrix_block_mean( dist_matrix, x_base, y_base, x_len, y_len, -1 );
            if(value != -1)
            {
                if(value < dist_1)
                {
                    if(color_1 != "NULL")
                    {
                        heatmap_matrix.at(idx).at(idy).fill = true;
                        heatmap_matrix.at(idx).at(idy).fill_color = QString::fromStdString(color_1);
                    }
                }
                else if(value < dist_2)
                {
                    if(color_2 != "NULL")
                    {
                        heatmap_matrix.at(idx).at(idy).fill = true;
                        heatmap_matrix.at(idx).at(idy).fill_color = QString::fromStdString(color_2);
                    }
                }
                else if(value < dist_3)
                {
                    if(color_3 != "NULL")
                    {
                        heatmap_matrix.at(idx).at(idy).fill = true;
                        heatmap_matrix.at(idx).at(idy).fill_color = QString::fromStdString(color_3);
                    }
                }else{
                    if(color_4 != "NULL")
                    {
                        heatmap_matrix.at(idx).at(idy).fill = true;
                        heatmap_matrix.at(idx).at(idy).fill_color = QString::fromStdString(color_4);
                    }
                }
            }
            else{
                heatmap_matrix.at(idx).at(idy).fill = true;
                heatmap_matrix.at(idx).at(idy).fill_color = Qt::black;
            }
        }
}




























