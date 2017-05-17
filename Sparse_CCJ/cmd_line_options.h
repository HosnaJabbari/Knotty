#pragma once

#include <string>
#include <fstream>
#include <iostream>

typedef struct cmd_line_args
{
private:
    bool changeable = true;

    bool use_sparse_ = true; // use sparse pseudo loop
    bool use_garbage_collection_ = true;
    char print_trace_arrow_info_ = 0;
    char print_candidate_list_info_ = 0;

public:
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

extern cmd_line_args cmd_line_options;
//extern shape_info shape;


