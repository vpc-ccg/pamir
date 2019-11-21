#ifndef __LOGGER__
#define __LOGGER__

#include <string>
#include <cstdarg>
#include <memory>
#include <cstdio>

using namespace std;

void log_init(string);

void log_close();

void log(const char *, ...);

//template<size_t BufferSize>
class Logger {
private:
    static const size_t BufferSize = 5 * 1024 * 1024;
    FILE *debugf;
    FILE *infof;
    FILE *errorf;
    char info_buf[2 * BufferSize] = {0};
    char error_buf[2 * BufferSize] = {0};
    size_t info_buf_size = 0;
    size_t error_buf_size = 0;

    /**
     * Flushes the current log buffer
     * @param fout file handler to flush
     * @param buffer_size current size of the used buffer
     * @param buffer current buffer
     */
    static void flush(FILE *fout, size_t &buffer_size, char *buffer) {
        fwrite(buffer, sizeof(char), buffer_size, fout);
        fflush(fout);
        buffer_size = 0;
    }

    /**
     * Flushes the buffer if it is full
     * @param fout fout file handler to dump
     * @param buffer_size current size of the used buffer
     * @param buffer current buffer
     */
    static void dump(FILE *fout, size_t &buffer_size, char *buffer) {
        if (buffer_size >= BufferSize) {
            flush(fout, buffer_size, buffer);
        }
    }

    Logger(const std::string &error_str, const std::string &info_str) : debugf(stderr) {
        errorf = fopen(error_str.c_str(), "w+");
        if (errorf == NULL) {
            fputs("Could not the open error stream\n", stderr);
            exit(-1);
        }
        infof = fopen(info_str.c_str(), "w+");

        if (infof == NULL) {
            fputs("Could not open the info stream\n", stderr);
            exit(-1);
        }
    }

    Logger(const std::string &error_str, const std::string &info_str, const std::string &debug_str) : Logger(error_str,
                                                                                                             info_str) {
        debugf = fopen(debug_str.c_str(), "w+");
        if (debugf == NULL) {
            fputs("Could not open the debug stream\n", stderr);
            exit(-1);
        }
    }

    Logger() : debugf(stderr), infof(stdout), errorf(stderr) {}

public:
    ~Logger() {
        flush(errorf, error_buf_size, error_buf);
        flush(infof, info_buf_size, info_buf);

        if (debugf != stderr) {
            fclose(debugf);
        }
        if (infof != stdout) {
            fclose(infof);
        }
        if (errorf != stderr) {
            fclose(errorf);
        }
    }


    Logger(const Logger &) = delete;

    Logger &operator=(const Logger &) = delete;

    Logger(Logger &&) = delete;

    Logger &operator=(Logger &&) = delete;

    /**
     * Get an instance of the logger object.
     * @param error Filename for error stream or "" for stderr
     * @param info Filename for info stream or "" for stdout
     * @param debug Filename for debug stream or "" for stderr
     * @return Reference to Logger object instance.

    static auto &instance(const std::string &error, const std::string &info, const std::string &debug) {
        static const std::unique_ptr <Logger> logger{new Logger<BufferSize>(error, info, debug)};
        return *logger;
    }
    */
    /**
     * Get an instance of the logger object using default stderr and stdout streams
     * @return Reference to Logger object instance.
     */
    static auto &instance() {
        static const std::unique_ptr <Logger> logger{new Logger()};
        //static Logger<BufferSize> logger;
        return *logger;
    }

    Logger &set_info(const std::string &str){
        infof = fopen(str.c_str(), "w+");
        if (infof == NULL) {
            fputs("Could not open the info stream\n", stderr);
            exit(-1);
        }
        return *this;
    }

    Logger &set_debug(const std::string &str){
        debugf = fopen(str.c_str(), "w+");
        if (debugf == NULL) {
            fputs("Could not open the debug stream\n", stderr);
            exit(-1);
        }
        return *this;
    }


    Logger &set_error(const std::string &str){
        errorf = fopen(str.c_str(), "w+");
        if (errorf == NULL) {
            fputs("Could not open the error stream\n", stderr);
            exit(-1);
        }
        return *this;
    }

    /**
     *  Logger command to log the debug level information
     * @param format
     * @param ...
     */
    void debug(const char *format, ...) {
#ifdef DEBUG
        va_list args;
        va_start(args, format);
        fprintf(debugf, format, args);
        fflush(debugf);
        va_end(args);
#endif
    }

    /**
     * Logger command to log the info level information
     * @param format
     * @param ...
     */
    void info(const char *format, ...) {
        va_list args;
        va_start(args, format);
        info_buf_size += vsprintf(info_buf + info_buf_size, format, args);
        va_end(args);
        dump(infof, info_buf_size, info_buf);
    }

    /**
     * Logger command to log the error level information
     * @param format
     * @param ...
     */
    void error(const char *format, ...) {
        va_list args;
        va_start(args, format);
        error_buf_size += vsprintf(error_buf + error_buf_size, format, args);
        va_end(args);
        dump(errorf, error_buf_size, error_buf);
    }

};

#endif
