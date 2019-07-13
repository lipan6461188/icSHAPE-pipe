
#ifndef PARIS_PLOT_H
#define PARIS_PLOT_H

#ifndef NO_QT

#include <QPainter>
#include <QPdfWriter>
#include <QDebug>

#endif

#include <cmath>
#include <algorithm>
#include <iostream>

#include "shape.h"
#include "pan_type.h"
#include "exceptions.h"
#include "param.h"
#include "align.h"

namespace pan{

using std::string;
using std::vector;
using std::unordered_map;
using std::istream;
using std::ostream;
using std::istringstream;
using std::pair;
using std::ifstream;
using std::runtime_error;
using std::min;

enum ColorMapType
{
    ENERGY
};

#ifndef NO_QT

struct Block
{
    double value; // raw value

    bool border=false;
    bool fill=false;

    QColor fill_color;
    QColor border_color;

    double border_width;
};

struct Diagonal_Rect
{
    double start = 0;
    double end = 0;

    QColor fill_color;
    QColor border_color;

    double border_width;
};

struct Heatmap_Matrix
{
    using size_type = Matrix<Block>::size_type;

    Matrix<Block> matrix;
    //vector<Diagonal_Rect> diagonal_rect;

    size_type size()const { return matrix.size(); }
};

//using Heatmap_Matrix = Matrix<Block>;

/*
void plot_heatmap( const Heatmap_Matrix &matrix, 
                    const string &out_pdf_file, 
                    double factor,
                    uLONG start_1=1, uLONG end_1=-1UL,
                    uLONG start_2=1, uLONG end_2=-1UL);
*/

void plot_intergrated_heatmap(const Heatmap_Matrix &matrix, 
            const string &out_pdf_file,
            const double factor,
            uLONG start_1=1, uLONG end_1=-1UL,
            uLONG start_2=1, uLONG end_2=-1UL,
            const string genome="", const Chr_Shape shape=Chr_Shape());

void plot_heatmap(  QPainter *p_painter,
                    const qreal width, const qreal height,
                    const Heatmap_Matrix &matrix, 
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2);

void plot_heatmap_sticks(  QPainter *p_painter,
                    const qreal width, const qreal height,
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2,
                    const bool up=true, const bool left=true, 
                    const bool down=true, const bool right=true);

void plot_labels(QPainter *p_painter,
                const qreal width, const qreal height,
                const QString &label_ver, const QString &label_her, 
                double bias, 
                const bool up=true, const bool left=true,
                const bool down=true, const bool right=true);


void plot_heatmap_sequence( QPainter *p_painter,
                            const qreal width, const qreal height,
                            const uLONG start_1, const uLONG end_1,
                            const uLONG start_2, const uLONG end_2,
                            const string &seq_1, const string &seq_2);

void plot_heatmap_shape_circle( QPainter *p_painter,
                            const qreal width, const qreal height,
                            const uLONG start_1, const uLONG end_1,
                            const uLONG start_2, const uLONG end_2,
                            const Chr_Shape &shape);

template<typename T>
void matrix_to_heatmap(     const Matrix<T> &matrix,
                            Heatmap_Matrix &heatmap_matrix,
                            const T min_v,
                            const T max_v,
                            const string &color="#C44E52")
{
    using size_type = typename Matrix<T>::size_type;

    init_matrix(heatmap_matrix.matrix, matrix.size());

    for(size_type idx=0; idx<matrix.size(); idx++)
        for(size_type idy=0; idy<matrix.size(); idy++)
        {
            auto value = matrix.at(idx).at(idy);
            if(value >= min_v)
            {
                heatmap_matrix.matrix.at(idx).at(idy).fill = true;
                heatmap_matrix.matrix.at(idx).at(idy).value = value;

                QColor cur_color( QString::fromStdString(color) );
                int alpha = 0;
                
                if(value >= max_v)
                    alpha = 255;
                else
                    alpha = min(int(1.0 * (value - min_v) / (max_v - min_v) * 255 + 0), 255);

                cur_color.setAlpha(alpha);
                heatmap_matrix.matrix.at(idx).at(idy).fill_color = cur_color;
            }
        }
}


void combine_heatmap(   const Heatmap_Matrix &heatmap_matrix_1,
                        const Heatmap_Matrix &heatmap_matrix_2,
                        Heatmap_Matrix &combine_heatmap_matrix );



void dark_gap_alignment_region( Heatmap_Matrix &heatmap,
                                const Multi_Align &alignment, 
                                const string &chromosome_name,
                                const string &color="#D1CFCF" );

#endif


#ifndef NO_STDOUT_PLOT

namespace STD_Color
{
    static Color::Modifier RED(Color::FG_RED);
    static Color::Modifier GREEN(Color::FG_GREEN);
    static Color::Modifier BLUE(Color::FG_BLUE);
    static Color::Modifier MAGENTA(Color::FG_MAGENTA);
    static Color::Modifier WHITE(Color::FG_WHITE);
    static Color::Modifier DEF(Color::FG_DEFAULT);
};

struct Color_Control
{
    double level_up_1 = 10;
    double level_up_2 = 50;
    double level_up_3 = 100;

    double level_dn_1 = 10;
    double level_dn_2 = 50;
    double level_dn_3 = 100;

    const Color::Modifier color_1 = STD_Color::BLUE; // level_1 - level_2
    const Color::Modifier color_2 = STD_Color::GREEN; // level_2 - level_3
    const Color::Modifier color_3 = STD_Color::RED; // > level_3

    bool valid() const;
    void print_legend(const string block_unicode = "\u25A0") const;
};

void print_heatmap( const Matrix<double> &matrix, 
                    const Color_Control &cc,
                    bool white_diagonal = false,
                    const uLONG indent = 0,
                    const string block_unicode = "\u25A0",
                    const Color::Modifier domain_color = STD_Color::WHITE,
                    const RegionArray &upper_region = RegionArray(),
                    const RegionArray &lower_region = RegionArray(),
                    const char sdomain_char = '*',
                    const bool domain_fill = false);

#endif




}
#endif