#ifndef __SEQT_HPP__
#define __SEQT_HPP__

#include <boost/compute.hpp>

#include <map>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>

using namespace boost::compute;

using std::map;
using std::ostream, std::wostream, std::cout, std::cerr, std::endl;
using std::ostream_iterator;


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
    kernel _find_nexts_kernel;
    kernel _collect_finds_kernel;
    kernel _mark_new_or_increment_count_kernel;
    kernel _initialize_newly_found_sequences_kernel;
    kernel _make_pair_constant_second_kernel;

    vector<long> _counts; // how many times have we seen this?

    vector<long> _lengths; // for sequences
    vector<long> _nexts; // for sequences
    vector<long> _prevs; // for sequences

    // vector<long> _depths; // for sets;
    // vector<long> _lefts; // for sets
    // vector<long> _rights; // for sets

    // 0 means not tracked
    // -1 means tracked at one character previous, -2 is two, etc.
    vector<long> _tracked;

    long _total;  

    map<wchar_t, long> _char_index;
    map<long, wchar_t> _index_char;

    long get_char_index(wchar_t c);

    vector<long> pack(vector<long> & data, function<long(long)> pred);
    void find_nexts(vector<long> & sorted_lengths, vector<long> & nexts_begin, vector<long> & nexts_end);
    void collect_finds(vector<long> & nexts_begin, vector<long> & scratch, vector<long> & current, vector<long2_> & found);
    void mark_new_or_increment_count(vector<long2_> & found, vector<long> & scratch);
    void initialize_newly_found_sequences(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value);
    vector<long2_> make_pair_constant_second(vector<long> & first, long second);

    long increment_and_add_new_finds(vector<long2_> & found, long tracked_value);

    template<typename T>
    void print(std::wostream & os, vector<T> & v); 

    void read(wchar_t c);
    void print_all(std::wostream & os);

    void print_by_index(std::wostream & os, long i, std::vector<long> const & prev, std::vector<long> const & next,
        std::map<long, std::wstring> & cache);
    std::wstring write_by_index(long i, std::vector<long> const & prev, std::vector<long> const & next,
        std::map<long, std::wstring> & cache);
    seqt();
}; 

    


#endif // __SEQT_HPP__