
#ifndef PARAM_H
#define PARAM_H

#include "pan_type.h"
#include <ctime>

namespace pan
{

void push_params(int &argc, char *argv[], const string &param);
std::string param_string(const int argc, char *argv[]);

#include <ostream>
namespace Color {
    enum Code {
    	FG_BLACK 			= 30,
        FG_RED      		= 31,
        FG_GREEN    		= 32,
        FG_YELLOW 			= 33,
        FG_BLUE     		= 34,
  		FG_MAGENTA 			= 35, 
		FG_CYAN 			= 36, 
		FG_LIGHT_GRAY 		= 37, 
        FG_DEFAULT  		= 39,
        BG_RED      		= 41,
        BG_GREEN    		= 42,
        BG_BLUE     		= 44,
        BG_DEFAULT  		= 49,		
		FG_DARK_GRAY 		= 90, 
		FG_LIGHT_RED 		= 91, 
		FG_LIGHT_GREEN 		= 92, 
		FG_LIGHT_YELLOW 	= 93, 
		FG_LIGHT_BLUE 		= 94, 
		FG_LIGHT_MAGENTA 	= 95, 
		FG_LIGHT_CYAN 		= 96, 
		FG_WHITE 			= 97
    };
    
    class Modifier {
        Code code;
    public:
        Modifier(Code pCode) : code(pCode) {}
        friend std::ostream&
        operator<<(std::ostream& os, const Modifier& mod) {
            return os << "\033[" << mod.code << "m";
        }
    };
}

namespace ColorArray
{
    // https://ws4.sinaimg.cn/large/006tNc79ly1fnlw1cnjlyj30dn02ndfo.jpg
    const StringArray Deep = { "#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#64b5cd" };

    // https://ws2.sinaimg.cn/large/006tNc79ly1fnlw1ti854j30dh02la9w.jpg
    const StringArray Muted = { "#4878cf", "#6acc65", "#d65f5f", "#b47cc7", "#c4ad66", "#77bedb" };

    // https://ws1.sinaimg.cn/large/006tNc79ly1fnlw2aeb8jj30dk02ja9w.jpg
    const StringArray Pastel = { "#92c6ff", "#97f0aa", "#ff9f9a", "#d0bbff", "#fffea3", "#b0e0e6" };

    // https://ws1.sinaimg.cn/large/006tNc79ly1fnlw2lkhpoj30dj02rdfo.jpg
    const StringArray Bright = { "#003fff", "#03ed3a", "#e8000b", "#8a2be2", "#ffc400", "#00d7ff" };

    // https://ws3.sinaimg.cn/large/006tNc79ly1fnlw2xddkfj30dc02pa9w.jpg
    const StringArray Dark = { "#001c7f", "#017517", "#8c0900", "#7600a1", "#b8860b", "#006374" };

    // https://ws3.sinaimg.cn/large/006tNc79ly1fnlw39phijj30di02oa9w.jpg
    const StringArray ColorBlind = { "#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56b4e9" };
}


namespace FILE_FORMAT
{
    enum file_format_type
    {
        SAM_FILE,
        BAM_FILE,
        FASTA_FILE,
        FASTQ_FILE,
        STOCKHOLM_FILE,
        CM_FILE,

        // I defined type
        MATRIX_FILE,
        DG_FILE,

        // UNKNOWM
        UNKNOWM
    };

    string fft_to_string(const file_format_type &fft);

    file_format_type guess_file_type(const string &file_name);
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime();




}
#endif