#pragma once

#include <string>
#include <fstream>
#include <iostream>

typedef struct cmd_line_args
{
private:
    int seq_len_;

    bool changeable = true;

    bool use_sparse_ = true; // use sparse pseudo loop
    bool use_garbage_collection_ = true;
    bool use_shape_data_ = false;
    std::string shape_file_;
    char print_trace_arrow_info_ = 0;
    char print_candidate_list_info_ = 0;

public:
    // prepares filename to be stored in shape_file
    void set_shape_file(std::string filename);

    void set_sequence_length(int length) {
        if (changeable)
            seq_len_ = length;
    }

    void set_use_sparse(bool use_sparse) {
        if (changeable)
            use_sparse_ = use_sparse;
    }

    void set_use_garbage_collection(bool use_garbage_collection) {
        if (changeable)
            use_garbage_collection_ = use_garbage_collection;
    }

    void set_print_trace_arrow_info(char ptai) {
        if (changeable)
            print_trace_arrow_info_ = ptai;
    }

    void set_print_candidate_list_info(char pcli) {
        if (changeable)
            print_candidate_list_info_ = pcli;
    }

    void set_done() {
        changeable = false;
    }

    const int sequence_length() const {
        return seq_len_;
    }

    const std::string shape_file() const{
        return shape_file_;
    }

    const bool use_shape_data() const {
        return use_shape_data_;
    }

    const bool use_sparse() const {
        return use_sparse_;
    }

    const bool use_garbage_collection() const {
        return use_garbage_collection_;
    }

    const char print_trace_arrow_info() const {
        return print_trace_arrow_info_;
    }

    const char print_candidate_list_info() const {
        return print_candidate_list_info_;
    }

}cmd_line_args;

class shape_info
{
private:
    /// TODO should maybe be changed/changeable
    /// Also, they should be with the other values, not here
    float b_ = -0.6; // intercept
    float m_ = 1.8; // slope
    double *data_;

public:
    shape_info() {}

    // based on http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
    bool is_number(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end()
                && (std::isdigit(*it) || *it == '.' || *it == '-' || *it == '+')
              ) ++it;
        return !s.empty() && it == s.end();
    }

    void set_data(std::string filename);
    void set_b(float b) { b_ = b; }
    void set_m(float m) { m_ = m; }

    double data(int index) { return data_[index]; }
    float b() { return b_; }
    float m() { return m_; }

};

extern cmd_line_args cmd_line_options;
extern shape_info shape;


