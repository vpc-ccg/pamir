#ifndef __LOGGER__
#define __LOGGER__

#include<string>

using namespace std;

void log_init(string);

void log_close();

void log(const char *, ...);

//Class Logger{
//private:
//    FILE *debug;
//    FILE *info;
//    FILE *error;
//    char *info_buf;
//    char *error_buf;
//    int info_buf_size;
//    int error_buf_size;
//
//    void dump (FILE *fout, int *buffer_size, char *buffer){
//        if (*buffer_size > 5 * 1024 * 1024) {
//            fwrite(buffer, sizeof(char), *buffer_size, output);
//            fflush(output);
//            *buffer_size = 0;
//        }
//    }
//    Logger ()
//
//public:
//    void static getLogger ()
//    void static getLogger (string, string, string)
//    /**
//     *  Logger command to log the debug level information
//     * @param format
//     * @param ...
//     */
//    void debug(const char *format, ...) {
//#ifdef DEBUG
//        va_list args;
//        va_start(args, format);
//        fprintf(debug, format, args);
//        fflush(debug);
//        va_end(args);
//#endif
//    }
//
//    /**
//     * Logger command to log the info level information
//     * @param format
//     * @param ...
//     */
//    void info(const char *format, ...) {
//        va_list args;
//        va_start(args, format);
//        info_buf_size += vsprintf(info_buf + info_buf_size, format, args);
//        va_end(args);
//        dump (info, debug_buf_size, debug_buf);
//    }
//
//    /**
//     * Logger command to log the error level information
//     * @param format
//     * @param ...
//     */
//    void error(const char *format, ...) {
//        va_list args;
//        va_start(args, format);
//        error_buf_size += vsprintf(error_buf + error_buf_size, format, args);
//        va_end(args);
//        dump (error, error_buf_size, error_buf);
//    }
//
//};
//
//
//
#endif
