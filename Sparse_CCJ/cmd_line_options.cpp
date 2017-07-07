#include "cmd_line_options.h"
#include <stdlib.h>    
#include <string>

cmd_line_args cmd_line_options;
shape_info shape = shape_info();

// prepares filename to be stored in shape_file
void cmd_line_args::set_shape_file(std::string filename) {
    if (changeable) {
        // cut off first 7 characters ("-shape=")
        filename = filename.substr(7,filename.length()-6);

        auto front = filename.front();
        auto back = filename.back();

        // cut quote from beginning and end if they exist
        if (front == '\"')
            filename.erase(front);
        if (back == '\"')
            filename.erase(back);

        shape_file_ = filename;
        use_shape_data_ = true;

        shape.set_data(shape_file_);
    }
}

void shape_info::set_data(std::string filename) {
        // open the file
        std::ifstream infile;
        infile.open(filename, std::ifstream::in);

        // size of data is sequence length
        int seq_len = cmd_line_options.sequence_length();
        data_ = new double[seq_len+1]; // data_[0] is not used - starts from 1

        // first line may be the sequence
        std::string input;
        char first_char;

        // continue until there is a line that is a number
        while (infile >> input && !is_number(input)) {}

        if (!is_number(input)) {
            printf("SHAPE data file error: file has no numbers in it\n");
            exit(-1);
        }

        // the first number is in input
        // it starts from 1
        data_[1] = stof(input);
        // go through where each word is a shape data number
        int i = 2;
        while (infile >> input) {
            if (i > seq_len) {
                printf("SHAPE data file error: length greater than sequence length (%d compared to %d)\n",i,seq_len);
                exit(-1);
            }
            if (!is_number(input)) {
                printf("SHAPE data file error: line after start is not a number\n",i,seq_len);
                exit(-1);
            }

            data_[i] = stof(input);
            ++i;
        }

        if (i-1 != seq_len) {
            printf("SHAPE data file error: length less than sequence length (%d compared to %d)\n",i-1,seq_len);
            exit(-1);
        }

        infile.close();
    }
