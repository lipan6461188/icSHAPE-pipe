
#include "paris_plot.h"
#include <exception>

using namespace std;

namespace pan{

#ifndef NO_QT

void plot_intergrated_heatmap(
            const Heatmap_Matrix &matrix,
            const string &out_pdf_file,
            const double factor,
            uLONG start_1, uLONG end_1,
            uLONG start_2, uLONG end_2,
            const string genome, const Chr_Shape shape)
{
    uLONG matrix_size = matrix.matrix.size();

    /* Check Range */
    if(start_1 > matrix_size or start_1 == 0)
        throw runtime_error("Out of Matrix Range");

    if(end_1 >= matrix_size)
        end_1 = matrix_size;

    if(start_1 > end_1)
        throw runtime_error("Bad Start End Range");


    if(start_2 > matrix_size or start_2 == 0)
        throw runtime_error("Out of Matrix Range");

    if(end_2 >= matrix_size)
        end_2 = matrix_size;

    if(start_2 > end_2)
        throw runtime_error("Bad Start End Range");

    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    double suggest_width = plot_matrix_width * factor * 5;
    double suggest_height = plot_matrix_height * factor * 5;
    
    QPdfWriter pdf_writer(QString::fromStdString(out_pdf_file));
    //pdf_writer.newPage();
    pdf_writer.newPage();
    pdf_writer.setResolution(1200);
    pdf_writer.setPageSize( QPageSize( QSizeF(suggest_width, suggest_height), QPageSize::Millimeter, "", QPageSize::ExactMatch ) );
    pdf_writer.setPageMargins(QMarginsF(100,100,100,100));
    
    QPainter painter( &pdf_writer );
    painter.setRenderHint( QPainter::Antialiasing );

    double width = painter.device()->width();
    double height = painter.device()->height();
    //cout << width << "\t" << height << endl;

    plot_heatmap( &painter, width, height, matrix, start_1, end_1, start_2, end_2);
    if(genome.empty() and shape.chr_id.empty())
        plot_heatmap_sticks( &painter, width, height, start_1, end_1, start_2, end_2, true, true, true, true);
    else
        plot_heatmap_sticks( &painter, width, height, start_1, end_1, start_2, end_2, false, false, true, true);


    if(not shape.chr_id.empty())
    {
        clog << "plot shape..." << endl;
        plot_heatmap_shape_circle( &painter, width, height, start_1, end_1, start_2, end_2, shape);
    }

    if(genome.size() != 0)
    {
        plot_heatmap_sequence( &painter, width, height, start_1, end_1, start_2, end_2, 
                    genome.substr(start_1-1, end_1-start_1+1), genome.substr(start_2-1, end_2-start_2+1));
    }


    /*

    double x_step = width / plot_matrix_width;
    double y_step = height / plot_matrix_height;

    double short_stick_len = (x_step+y_step)/2;
    double long_stick_len = 2*short_stick_len;
    */

    //double bias = QFontMetrics(QFont("Helvetica Neue")).height() + long_stick_len;

    //plot_labels(&painter, width, height, "Ver", "Hev", 2*long_stick_len, true, true, true, true);

}


void plot_heatmap(  QPainter *p_painter,
                    const qreal width, const qreal height,
                    const Heatmap_Matrix &matrix, 
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2)
/* start-end are 1-based and closed */
{
    using size_type = Matrix<Block>::size_type;

    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;


    double x_step = width / plot_matrix_width;
    double y_step = height / plot_matrix_height;

    QBrush fill_brush(Qt::SolidPattern);
    QPen border_pen(Qt::SolidLine);

    p_painter->setPen(Qt::NoPen);
    for(size_type x_idx=start_1-1;x_idx<end_1;x_idx++)
        for(size_type y_idx=start_2-1;y_idx<end_2;y_idx++)
        {
            const auto &cur_block = matrix.matrix.at(x_idx).at(y_idx);

            double x_start = x_step * (x_idx-start_1+1);
            double y_start = y_step * (y_idx-start_2+1);
            //QRectF rect(y_start, x_start, y_step, x_step);
            QRectF rect(x_start, y_start, x_step, y_step);

            if(cur_block.fill)
            {
                fill_brush.setColor(cur_block.fill_color);
                p_painter->setBrush(fill_brush);
                p_painter->drawRect( rect );
            }
        }

    p_painter->setBrush(Qt::NoBrush);
    for(size_type x_idx=start_1-1;x_idx<end_1;x_idx++)
        for(size_type y_idx=start_2-1;y_idx<end_2;y_idx++)
        {
            const auto &cur_block = matrix.matrix.at(x_idx).at(y_idx);

            double x_start = x_step * (x_idx-start_1+1);
            double y_start = y_step * (y_idx-start_2+1);
            //QRectF rect(y_start, x_start, y_step, x_step);
            QRectF rect(x_start, y_start, x_step, y_step);

            if(cur_block.border)
            {
                border_pen.setWidthF(cur_block.border_width);
                border_pen.setColor(cur_block.border_color);
                p_painter->setPen(border_pen);
                p_painter->drawRect( rect );
            }
        }

}

void plot_heatmap_sticks(  QPainter *p_painter,
                    const qreal width, const qreal height,
                    const uLONG start_1, const uLONG end_1,
                    const uLONG start_2, const uLONG end_2,
                    const bool up, const bool left, 
                    const bool down, const bool right)
{

    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    //double suggest_width = plot_matrix_width * factor;
    //double suggest_height = plot_matrix_height * factor;

    double x_step = width / plot_matrix_width;
    double y_step = height / plot_matrix_height;

    QPen stick_pen;

    stick_pen.setColor(Qt::black);
    stick_pen.setWidthF(1.0);
    p_painter->setPen(stick_pen);

    QFont label_font("Helvetica Neue");
    p_painter->setFont(label_font);

    double short_stick_len = (x_step+y_step)/2;
    double long_stick_len = 2*short_stick_len;


    for(uLONG x=0; x<=plot_matrix_height; x++)
    {
        double left_x_start = 0;
        double y_start = y_step * x;

        double right_x_start = width;
        
        uLONG stick_label = x + start_2 - 1;

        if(stick_label % 10 == 0)
        {
            if(left)
            {
                p_painter->drawLine(left_x_start, y_start, left_x_start-long_stick_len, y_start);
                p_painter->drawText(QRectF(left_x_start-long_stick_len-1500, y_start-200, 1500, 400), Qt::AlignRight | Qt::AlignVCenter, QString::number(stick_label));
            }
            if(right)
            {
                p_painter->drawLine(right_x_start, y_start, right_x_start+long_stick_len, y_start);
                p_painter->drawText(QRectF(right_x_start+long_stick_len, y_start-200, 1500, 400), Qt::AlignLeft | Qt::AlignVCenter, QString::number(stick_label));
            }
        }
        else
        {
            if(left)
                p_painter->drawLine(left_x_start, y_start, left_x_start-short_stick_len, y_start);
            if(right)
                p_painter->drawLine(right_x_start, y_start, right_x_start+short_stick_len, y_start);
        }
    }

    for(uLONG x=0; x<=plot_matrix_width; x++)
    {
        double top_y_start = 0;
        double x_start = x_step * x;

        double bot_y_start = height;

        uLONG stick_label = x + start_1 - 1;

        if(stick_label % 10 == 0)
        {
            if(up)
            {
                p_painter->drawLine(x_start, top_y_start, x_start, top_y_start-long_stick_len);
                p_painter->drawText(QRectF(x_start-500, top_y_start-long_stick_len-500, 1000, 500), Qt::AlignBottom | Qt::AlignHCenter, QString::number(stick_label));
            }

            if(down)
            {
                p_painter->drawLine(x_start, bot_y_start, x_start, bot_y_start+long_stick_len);
                p_painter->drawText(QRectF(x_start-500, bot_y_start+long_stick_len, 1000, 500), Qt::AlignTop | Qt::AlignHCenter, QString::number(stick_label));
            }
        }
        else
        {
            if(up)
                p_painter->drawLine(x_start, top_y_start, x_start, top_y_start-short_stick_len);
            if(down)
                p_painter->drawLine(x_start, bot_y_start, x_start, bot_y_start+short_stick_len);
        }

    }
}


void plot_labels(QPainter *p_painter,
                const qreal width, const qreal height,
                const QString &label_ver, const QString &label_her, 
                double bias, 
                const bool up, const bool left,
                const bool down, const bool right)
{
    QPen label_pen;
    QFont label_font("Helvetica Neue");

    label_pen.setColor(Qt::black);
    p_painter->setPen(label_pen);
    p_painter->setFont(label_font);

    double x_center = width / 2;
    double y_center = height / 2;

    double y_up = 0 - bias;
    double y_down = height + bias;
    //double x_left = 0 - long_stick_len * 2;
    double x_right = width + bias;

    if(up)
        p_painter->drawText(x_center-500, y_up-500, 1000, 500, Qt::AlignBottom | Qt::AlignHCenter, label_ver);

    if(down)
        p_painter->drawText(x_center-500, y_down, 1000, 500, Qt::AlignTop | Qt::AlignHCenter, label_ver);

    if(left)
    {
        p_painter->rotate(-90);
        p_painter->drawText(-y_center-500, -bias-500, 1000, 500, Qt::AlignBottom | Qt::AlignHCenter, label_her);
        //p_painter->drawText(200, 100, 500, 1000, Qt::AlignLeft | Qt::AlignTop, "HAHAH");
        p_painter->rotate(90);
    }

    if(right)
    {
        p_painter->rotate(90);
        p_painter->drawText(y_center-500, 0-x_right-1000, 1000, 1000, Qt::AlignBottom | Qt::AlignHCenter, label_her);
        p_painter->rotate(-90);
    }
}


void plot_heatmap_sequence( QPainter *p_painter,
                            const qreal width, const qreal height,
                            const uLONG start_1, const uLONG end_1,
                            const uLONG start_2, const uLONG end_2,
                            const string &seq_1, const string &seq_2)
{
    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    double x_step = width / plot_matrix_width;
    double y_step = height / plot_matrix_height;

    QPen base_pen;
    base_pen.setColor(Qt::black);
    p_painter->setPen(base_pen);

    QFont base_font("Helvetica Neue");
    base_font.setPointSizeF(base_font.pointSizeF()-1);
    p_painter->setFont(base_font);


    for(uLONG x=0; x<plot_matrix_height; x++)
    {
        double left_x_start = 0;
        double y_start = y_step * x;

        QString base = QString::fromStdString( seq_2.substr(x,1) );
        p_painter->drawText(left_x_start-y_step, y_start, y_step, y_step, Qt::AlignHCenter | Qt::AlignVCenter, base);
    }

    for(uLONG x=0; x<plot_matrix_width; x++)
    {
        double top_y_start = 0;
        double x_start = x_step * x;

        QString base = QString::fromStdString( seq_1.substr(x,1) );
        //p_painter->drawText(x_start, top_y_start-500, x_step, 500, Qt::AlignHCenter | Qt::AlignBottom, base);
        p_painter->drawText(x_start, top_y_start-x_step, x_step, x_step, Qt::AlignHCenter | Qt::AlignVCenter, base);
    }
}

QColor color_map(qreal number_0_1, const ColorMapType type)
{
    switch(type)
    {
        case ENERGY:
            if( number_0_1==0.00 )
                return QColor("#4747B6");
            else if( number_0_1 < 0.20 )
                return QColor("#4747FF");
            else if( number_0_1 < 0.40 )
                return QColor("#1CFF47");
            else if( number_0_1 < 0.60 )
                return QColor("#FFFF47");
            else if( number_0_1 < 0.80 )
                return QColor("#FF4747");
            else
                return QColor("#B64747");
            break;
        default:
            return Qt::gray;
    }

}

void plot_heatmap_shape_circle( QPainter *p_painter,
                            const qreal width, const qreal height,
                            const uLONG start_1, const uLONG end_1,
                            const uLONG start_2, const uLONG end_2,
                            const Chr_Shape &shape)
{
    uLONG plot_matrix_width = end_1 - start_1 + 1;
    uLONG plot_matrix_height = end_2 - start_2 + 1;

    double x_step = width / plot_matrix_width;
    double y_step = height / plot_matrix_height;

    //p_painter->setPen(Qt::NoPen);

    for(uLONG x=0; x<plot_matrix_height; x++)
    {
        double left_x_start = 0;
        double y_start = y_step * x;

        float shape_value = shape.shape_score[start_2+x-1];
        if(shape_value == NULL_Score)
        {
            p_painter->setBrush(Qt::gray);
        }else{
            p_painter->setBrush( color_map(shape_value, ENERGY) );
        }

        p_painter->drawEllipse(left_x_start-y_step, y_start, y_step, y_step);
    }

    // up side
    for(uLONG x=0; x<plot_matrix_width; x++)
    {
        double top_y_start = 0;
        double x_start = x_step * x;

        float shape_value = shape.shape_score[start_1+x-1];
        if(shape_value == NULL_Score)
        {
            p_painter->setBrush(Qt::gray);
        }else{
            p_painter->setBrush( color_map(shape_value, ENERGY) );
        }

        p_painter->drawEllipse(x_start, top_y_start-x_step, x_step, x_step);
    }
}



void combine_heatmap(   const Heatmap_Matrix &heatmap_matrix_1,
                        const Heatmap_Matrix &heatmap_matrix_2,
                        Heatmap_Matrix &combine_heatmap_matrix )
{
    if(heatmap_matrix_1.size() != heatmap_matrix_2.size())
        throw Unexpected_Error("Different Size of 2 input ");

    Heatmap_Matrix::size_type matrix_size = heatmap_matrix_1.size();

    combine_heatmap_matrix = heatmap_matrix_1;
    
    for(uLONG idx=0; idx<matrix_size; idx++)
        for(uLONG idy=idx+1; idy<matrix_size; idy++)
            combine_heatmap_matrix.matrix[idy][idx] = heatmap_matrix_2.matrix[idx][idy]; //rm_2.at(idx).at(idy).count;
}


void dark_gap_alignment_region( Heatmap_Matrix &heatmap,
                                const Multi_Align &alignment, 
                                const string &chromosome_name,
                                const string &color )
{
    using namespace std::placeholders;
    auto matrix_size = heatmap.size();

    if(not alignment.has(chromosome_name))
        throw Unexpected_Error("Chromosome "+chromosome_name+" is not in alignment file");

    auto coor_covert = std::bind(&Multi_Align::raw_coor_to_align_coor, _1, chromosome_name, _2);

    uLONG align_length = alignment.length();
    uLONG raw_length = alignment.get_sto_record(chromosome_name).seq_length;
    double step_size = 1.0 * align_length / matrix_size;

    Matrix<uLONG> count_matrix;
    init_matrix(count_matrix, matrix_size);

    for(uLONG x=0; x<raw_length; x++)
        for(uLONG y=0; y<raw_length; y++)
        {
            uLONG x_convert = coor_covert(alignment, x) - 1;
            uLONG y_convert = coor_covert(alignment, y) - 1;
            x_convert /= step_size;
            y_convert /= step_size;
            count_matrix.at(x_convert).at(y_convert)++;
        }

    for(uLONG x=0; x<matrix_size; x++)
        for(uLONG y=0; y<matrix_size; y++)
        {
            if( count_matrix[x][y] < uLONG(0.5*step_size) )
            {
                heatmap.matrix[x][y].border = false;
                heatmap.matrix[x][y].fill = true;
                heatmap.matrix[x][y].fill_color = QColor( QString::fromStdString(color) );
            }    
        }
}

#endif

#ifndef NO_STDOUT_PLOT

bool Color_Control::valid() const
{
    if( level_dn_1 < level_dn_2 and level_dn_2 < level_dn_3 and level_up_1 < level_up_2 and level_up_2 < level_up_3 )
        return true;
    else
        return false;
}

void Color_Control::print_legend(const string block_unicode) const
{
    cout << "color legend: (right) "  << color_1 << block_unicode+" " << " " << level_up_1 << "~" << level_up_2 << "    "
                                    << color_2 << block_unicode+" " << " " << level_up_2 << "~" << level_up_3 << "    "
                                    << color_3 << block_unicode+" " << " >" << level_up_3  << STD_Color::DEF;
    cout << "     (left) "  << color_1 << block_unicode+" " << " " << level_dn_1 << "~" << level_dn_2 << "    "
                                    << color_2 << block_unicode+" " << " " << level_dn_2 << "~" << level_dn_3 << "    "
                                    << color_3 << block_unicode+" " << " >" << level_dn_3 << "\n" << STD_Color::DEF;
}

void print_heatmap( const Matrix<double> &matrix, 
                    const Color_Control &cc,
                    bool white_diagonal,
                    const uLONG indent,
                    const string block_unicode,
                    const Color::Modifier domain_color,
                    const RegionArray &upper_region,
                    const RegionArray &lower_region,
                    const char domain_char,
                    const bool domain_fill)
{
    auto matrix_size = matrix.size();
    const string indent_string(indent, ' ');

    if(not sorted_no_overlap(upper_region) or not sorted_no_overlap(lower_region))
        throw Unexpected_Error("unsorted or overlapped domain regions");
    if(not upper_region.empty() and upper_region.back().second > matrix_size)
        throw Unexpected_Error("domain region out of range");
    if(not lower_region.empty() and lower_region.back().second > matrix_size)
        throw Unexpected_Error("domain region out of range");

    uLONG up_first_domain_length = 0, dw_first_domain_length = 0;
    if(upper_region.size() != 0 and upper_region[0].first == 1)
        up_first_domain_length = upper_region[0].second;
    if(lower_region.size() != 0 and lower_region[0].first == 1)
        dw_first_domain_length = lower_region[0].second;

    cout << STD_Color::DEF << indent_string;
    // print head line
    //cout << "- ";
    for(uLONG idx=0; idx<=matrix_size; idx++)
    {
        cout << "-";
        if(idx <= up_first_domain_length)
            cout << domain_color << domain_char << STD_Color::DEF;
        else
            cout << " ";
    }

    cout << "\n";


    // print each line
    for(uLONG idx=0; idx<matrix_size; idx++)
    {
        cout << STD_Color::DEF << indent_string << "-";

        if(idx < dw_first_domain_length)
            cout << domain_color << domain_char << STD_Color::DEF;
        else
            cout << " ";

        for(uLONG idy=0; idy<matrix_size; idy++)
        {
            if(white_diagonal and idx == idy)
            {
                cout << STD_Color::WHITE << block_unicode;
                //continue;
            }else{
                if(idx < idy)
                // upper half
                {
                    if(matrix[idx][idy] < 0)
                    {
                        cout << STD_Color::MAGENTA << block_unicode;
                    }else if(matrix[idx][idy]<cc.level_up_1)
                    {
                        cout << " ";
                    }else if(matrix[idx][idy] < cc.level_up_2)
                    {
                        cout << cc.color_1 << block_unicode;//+" ";
                    }else if(matrix[idx][idy] < cc.level_up_3)
                    {
                        cout << cc.color_2 << block_unicode;//+" ";
                    }else
                    {
                        cout << cc.color_3 << block_unicode;//+" ";
                    }
                }else{
                // lower half
                    if(matrix[idx][idy] < 0)
                    {
                        cout << STD_Color::MAGENTA << block_unicode;
                    }else if(matrix[idx][idy]<cc.level_dn_1)
                    {
                        cout << " ";
                    }else if(matrix[idx][idy] < cc.level_dn_2)
                    {
                        cout << cc.color_1 << block_unicode;//+" ";
                    }else if(matrix[idx][idy] < cc.level_dn_3)
                    {
                        cout << cc.color_2 << block_unicode;//+" ";
                    }else if(matrix[idx][idy] < 0)
                    {
                        cout << STD_Color::MAGENTA << block_unicode;
                    }else
                    {
                        cout << cc.color_3 << block_unicode;//+" ";
                    }
                }
            }
            
            // print domain char            
            if(idx <= idy)
            // upper half
            {
                bool to_plot;
                if(domain_fill)
                    to_plot = std::any_of(upper_region.cbegin(), upper_region.cend(), [&](const Region &g){ return (g.first <= idx+1 and idy+1 <= g.second); });
                else
                    to_plot = std::any_of(upper_region.cbegin(), upper_region.cend(), [&](const Region &g){ return (g.first <= idx+1 and idy+1 == g.second) or (g.first == idx+1 and idy+1 <= g.second); });

                if(to_plot)
                    cout << domain_color << domain_char << STD_Color::DEF;
                else
                    cout << " ";
            }else if(idx > idy)
            // lower half
            {
                bool to_plot;
                if(domain_fill) 
                    to_plot = std::any_of(lower_region.cbegin(), lower_region.cend(), [&](const Region &g){ return (g.first <= idy+2 and idx+1 <= g.second); });
                else
                    to_plot = std::any_of(lower_region.cbegin(), lower_region.cend(), [&](const Region &g){ return (g.first <= idy+2 and idx+1 == g.second) or (g.first == idy+2 and idx+1 <= g.second); });

                if(to_plot)
                    cout << domain_color << domain_char << STD_Color::DEF;
                else
                    cout << " ";
            }

        }
        cout << STD_Color::DEF << "-\n";
    }
    cout << STD_Color::DEF;

    // print tail line
    cout << indent_string;
    for(uLONG idx=0; idx<matrix_size+1; idx++)
        cout << "- ";
    cout << "\n";
}

#endif




}