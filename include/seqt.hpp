#ifndef __SEQT_HPP__
#define __SEQT_HPP__

#include <boost/compute.hpp>

#include <map>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>


using std::map;
using std::ostream, std::wostream, std::cout, std::cerr, std::endl;
using std::ostream_iterator;

using boost::compute::device;
using boost::compute::context;
using boost::compute::command_queue;
using boost::compute::function;
using boost::compute::program;
using boost::compute::kernel;
using boost::compute::buffer_iterator;
using boost::compute::vector;
using boost::compute::long2_;

// BOOST_COMPUTE_FUNCTION(long, is_current, (long tracked), {
//     if(tracked == 0) 
//         return 1;
//     return 0;
// });

long calc_operational_size(long global_size, long local_size);

// ostream& operator<<(ostream & os, long2_ x) {
//     os << "(" << 
// }

class seqt { 
public:
    device _device;
    context _context;
    command_queue _queue;
    program _program;
    kernel _pack_kernel;
    kernel _scatter_value_kernel;
    kernel _find_nexts_kernel;
    kernel _collect_finds_kernel;
    kernel _mark_exists_kernel;
    kernel _initialize_newly_found_sequences_kernel;
    kernel _make_pair_constant_second_kernel;
    kernel _calculate_stats_kernel;
    kernel _is_sequence_significant_kernel;
    kernel _depends_on_sorted_list_kernel;
    kernel _depends_on_flagged;
    kernel _gather_seqs_kernel;

    vector<long> _counts; // how many times have we seen this?

    vector<long> _lengths; // for sequences
    vector<long2_> _seqs; // for sequences
    vector<long2_> _initial_seq_counts; // what were the counts of the children at the time this sequence was created?
    vector<long> _initial_characters_read; // how many characters were read at the time of this sequence creation?

    vector<float> _expected_counts;
    vector<float> _stddev_counts;
    vector<float> _significance;

    // vector<long> _depths; // for sets;
    // vector<long> _lefts; // for sets
    // vector<long> _rights; // for sets

    // 0 means not tracked
    // -1 means tracked at one character previous, -2 is two, etc.
    vector<long> _tracked;

    long _total;  
    long _characters_read;

    float _sigma;
    float _min_sigma; // what do we throw away during sleep?
    long _min_occurances;
    long _max_sequences_tracked;

    map<wchar_t, long> _char_index;
    map<long, wchar_t> _index_char;
    map<long, std::wstring> _seq_strings;

    long get_char_index(wchar_t c);

    vector<long> pack(vector<long> & data, function<long(long)> pred);
    void scatter_value(vector<long> & indices, long value, vector<long> & output);
    void find_nexts(vector<long> & sorted_lengths, vector<long2_> & nexts); 
    void collect_finds(vector<long2_> & nexts_begin, vector<long> & scratch, vector<long> & current, vector<long2_> & found);
    void mark_exists(vector<long2_> & found, vector<long> & scratch);
    void initialize_newly_found_sequences(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value);
    vector<long2_> find_nexts_by_length(vector<long> & current);
    vector<long2_> make_pair_constant_second(vector<long> & first, long second);
    void is_sequence_significant(vector<long2_> & seq, float sigma, long min_count, vector<long> & output);
    void depends_on_sorted_list(buffer_iterator<long> sorted_begin, buffer_iterator<long> sorted_end, 
                            buffer_iterator<long> begin, buffer_iterator<long> end,
                            vector<long> & output);
    void depends_on_flagged(vector<long> & flagged, vector<long> & output);
    void recreate_from_map(vector<long> & index);
    void calculate_stats();
    void gather_seqs(vector<long> & index, vector<long> & reverse_index, vector<long2_> & output);

    bool flag_existing(vector<long> & existing_indices, vector<long> & current_flag);
    long process_new_finds(vector<long2_> & finds, vector<long> & current_flag);

    long add_new_finds(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value);

    void remove_least_significant(long max_sequences);

    template<typename T>
    void print(std::wostream & os, vector<T> & v); 

    void read(wchar_t c);
    void print_all(std::wostream & os);

    void print_by_index(std::wostream & os, long i, std::vector<long2_> const & seqs,
        std::vector<long> const & counts, std::vector<float> & sig, std::map<long, std::wstring> & cache);

    std::wstring write_by_index(long i, std::vector<long2_> const & seqs,
        std::map<long, std::wstring> & cache);
        
    seqt();
}; 

    


#endif // __SEQT_HPP__