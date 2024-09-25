#include "seqt.hpp"

#include <cmath>
#include <iterator>
#include <algorithm>
#include <tuple>

#include "nth_element_struct.hpp"

#include <boost/compute/function.hpp>

using std::tuple;
using std::wcout;
using std::make_tuple;


using namespace boost::compute;


// number of standard deviations from the expected value ab is
float sequence_stddevs(long a, long b, long ab) {
    using std::sqrt;

    if(ab > a || ab > b || a + b == 0)
        return 0;

    float total = (float)(a + b);
    float pa = (float)a / total;
    float pb = (float)b / total;
    float pab = pa * pb;
    float expected_value = total * pab;
    float variation = total * pab * (1-pab);
    float stddev = sqrt(variation);

    return ((float)ab - expected_value) / stddev;
}

float sequence_pvalue(long a, long b, long ab) {
    using std::sqrt;
    using std::erfc;

    if(a == 0 && b == 0) {
        if(ab > 0) return 0.;
        return 0.5;
    }
    float total = (float)(a + b);
    float pa = (float)a / total;
    float pb = (float)b / total;
    float pab = pa * pb;
    float expected_value = total * pab;
    float variation = total * pab * (1-pab);

    float x = ((float)ab - expected_value) * sqrt(2 / variation);

    return 0.5 * erfc(x);
}


// struct significance_index {
//     double significance;
//     long index;

//     bool operator==(significance_index const & r) const {
//         return significance == r.significance && index == r.index;
//     }
// };

// std::wostream & operator<<(std::wostream & os, significance_index const & r) {
//     os << "{" << r.index << ":" << r.significance << "}";
//     return os;
// }

// BOOST_COMPUTE_ADAPT_STRUCT(significance_index, significance_index, (significance, index));

// BOOST_COMPUTE_FUNCTION(bool, less_significance_index, (significance_index a, significance_index b), {
//     if(a.significance < b.significance)
//         return true;
//     if(a.significance > b.significance)
//         return false;
//     return a.index < b.index;
// });

// BOOST_COMPUTE_FUNCTION(long, index_from_significance_index, (significance_index a), {
//     return a.index;
// });

// BOOST_COMPUTE_FUNCTION(significance_index, zip_significance_index, (float sig, long index), {
//     significance_index r;
//     r.significance = sig;
//     r.index = index;
//     return r;
// });

BOOST_COMPUTE_FUNCTION(long, flag_tracked, (long x), {
    if(x == 0)
        return 1;
    return x;
});

BOOST_COMPUTE_FUNCTION(long, unflag_tracked, (long x), {
    if(x == 1)
        return 0;
    return x;
});

BOOST_COMPUTE_FUNCTION(bool, long2_compare, (boost::compute::long2_ a, boost::compute::long2_ b), {
    if(a.x < b.x) return true;
    if(b.x < a.x) return false;
    if(a.y < b.y) return true;
    return false;
});

BOOST_COMPUTE_FUNCTION(long, pair_difference, (boost::compute::long2_ a), {
    return a.y - a.x;
});

BOOST_COMPUTE_FUNCTION(int, check_nan, (float a), {
    if(isnan(a)) 
        return 1;
    return 0;
});


long seqt::get_char_index(wchar_t c) {
    long index;

    auto i = _char_index.find(c);
    if(i == _char_index.end()) {
        index = _total;
        _char_index.insert({c, index});
        _index_char.insert({index, c});

        _total++;

        _counts.resize(_total, _queue);
        _initial_characters_read.resize(_total, _queue);
        _lengths.resize(_total, _queue);
        _seqs.resize(_total, _queue);
        _initial_seq_counts.resize(_total, _queue);
        _tracked.resize(_total, _queue);

        fill_n(_initial_characters_read.begin() + index, 1, _characters_read, _queue);
        fill_n(_lengths.begin() + index, 1, 1, _queue);

        // atoms have prev[i] and next[i] == 0
        fill_n(_seqs.begin() + index, 1, long2_(0,0), _queue);
        fill_n(_initial_seq_counts.begin() + index, 1, long2_(0,0), _queue);

        fill_n(_counts.begin() + index, 1, 0, _queue); // initialize count to zero because we haven't tracked this yet
        fill_n(_tracked.begin() + index, 1, 0, _queue);
    } else {
        index = i->second;// increment the count

        // don't update the count here, we'll do it in the main read loop
        // long c;
        // copy_n(_counts.begin() + index, 1, &c, _queue);
        // fill_n(_counts.begin() + index, 1, c + 1, _queue);
    }

    return index;
} 

void seqt::read(wchar_t c) {

    // std::wcout << "reading: " << c << endl;

    vector<long> current(0, _context);
    vector<long> scratch(0, _context);
    vector<long> existing_indices(0, _context);
    vector<long> new_find_indices(0, _context);
    vector<long2_> found(0, _context);
    vector<long2_> new_finds(0, _context);
    long found_count;
    long cnt;
    long existing_flagged;
    long new_count;
    long new_finds_added;

    // subtract one from all our tracked values
    transform(_tracked.begin(), _tracked.end(), _tracked.begin(), _1 - 1, _queue);

    // get the index of the current char
    // this will set the tracked value to 0 if it's new
    long index = get_char_index(c); 
    calculate_stats();
    _characters_read++;

    // flag our current sequences 
    vector<long> current_flag(_total, _context);
    transform(_tracked.begin(), _tracked.end(), current_flag.begin(), _1 == 0, _queue);
    
    // manually set the current flag to 1 for the index of
    // the current char
    fill_n(current_flag.begin() + index, 1, 1, _queue);

    bool do_remove_least_significant = false;

    for(;;) {
        // STEP 0: get all the current sequences
        current = pack(current_flag, _1 == 1);

        // create a list of pairs of potential sequences based on the length
        found = find_nexts_by_length(current);
        if(found.size() == 0)
            break; 

        // sort it to use it as an index
        sort(found.begin(), found.end(), long2_compare, _queue);

        // allocate memory and initialize to zero for marking new sequences as existing
        scratch = vector<long>(found.size(), _context);
        fill(scratch.begin(), scratch.end(), 0, _queue);

        // mark the indices in scratch for new ones with a 0, mark the found sequence for found ones
        mark_exists(found, scratch);

        // Collect all the new sequences and add them to our list to be tracked
        // now we have indices to all new sequences and existing ones
        new_find_indices = pack(scratch, _1 == 0);

        if(new_find_indices.size() + _total > _max_sequences_tracked) {
            do_remove_least_significant = true;
            new_finds_added = 0;
        } else {
            new_finds = vector<long2_>(new_find_indices.size(), _context);

            // gather all the new finds using the indices
            gather(new_find_indices.begin(), new_find_indices.end(), found.begin(), new_finds.begin(), _queue);
            new_finds_added = process_new_finds(new_finds, current_flag);
        }

        // first remove all the zeros, the remaining ones are the indices directly to our found sequences
        auto existing_end = remove(scratch.begin(), scratch.end(), 0, _queue);
        existing_indices = vector<long>(existing_end - scratch.begin(), _context);
        copy(scratch.begin(), existing_end, existing_indices.begin(), _queue);

        // now process our existing and new finds
        bool existing_were_flagged = flag_existing(existing_indices, current_flag);
        
        if(new_finds_added == 0 && !existing_were_flagged) {
            // we have new new indices and no existing ones flagged
            break;
        }
    }

    // increment the counts of all the flagged sequences
    // std::wcout << "current: ";
    // print(std::wcout, current);
    scratch = vector<long>(current.size(), _context);
    
    // increment counts for all the current sequences
    gather(current.begin(), current.end(), _counts.begin(), scratch.begin(), _queue);
    transform(scratch.begin(), scratch.end(), scratch.begin(), _1 + 1, _queue);
    scatter(scratch.begin(), scratch.end(), current.begin(), _counts.begin(), _queue);
    // std::wcout << "counts: ";
    // print(std::wcout, _counts);

    // mark all the current sequences as tracked=0
    scatter_value(current, 0, _tracked);

    if(do_remove_least_significant) {
        std::wcout << "\nremoving least significant" << endl;
        remove_least_significant(_total / 2);
    }
}

long seqt::process_new_finds(vector<long2_> & finds, vector<long> & current_flag) {
    if(finds.size() == 0)
        return 0;

    vector<long> scratch(finds.size(), _context);
    // filter the new sequences by significance
    // transform(_significance.begin(), _significance.end(), scratch2.begin(), scratch2.begin(), _2 == 0 && _1 > 5., _queue);
    is_sequence_significant(finds, _sigma, _min_occurances, scratch);

    // is it both significant and new?
    vector<long> new_find_indices = pack(scratch, _1 == 1);
    long new_count = new_find_indices.size();

    // std::wcout << "new_count after: " << new_count << endl;

    if(new_count > 0) {
        add_new_finds(new_find_indices, finds, 0);
        calculate_stats();
        
        long original_size = current_flag.size();
        current_flag.resize(_total, _queue);

        // all new sequences are flaged as current
        fill(current_flag.begin() + original_size, current_flag.end(), 1, _queue); 
    }

    return new_count;
}

bool seqt::flag_existing(vector<long> & existing_indices, vector<long> & current_flag) {
    if(existing_indices.size() == 0)
        return false;

    vector<long> scratch(existing_indices.size(), _context);

    // re-use scratch to see if these existing sequences have been flagged
    gather(existing_indices.begin(), existing_indices.end(), current_flag.begin(), scratch.begin(), _queue);

    long existing_flagged = existing_indices.size();
    reduce(scratch.begin(), scratch.begin() + existing_indices.size(), &existing_flagged, _queue);

    if(existing_flagged < existing_indices.size()) {
        // TODO:
        scatter_value(existing_indices, 1, current_flag);
        return true;
    }

    return false;
}

vector<long2_> seqt::find_nexts_by_length(vector<long> & current) {
    // STEP 1: sort all current sequences by length
    // allocate space for the current lengths

    vector<long2_> found(0, _context);

    vector<long> current_lengths(current.size(), _context);

    // gather the lengths of the current items into our allocated space
    gather(current.begin(), current.end(), _lengths.begin(), current_lengths.begin(), _queue);

    // sort the current indices by the lengths
    sort_by_key(current_lengths.begin(), current_lengths.end(), current.begin(), _queue);

    // STEP 2: Find potentially new sequences by looking how far back we saw
    //   a given sequence and marking in our sorted current lengths the possible
    //   subsequent sequences
    // now go through and find new sequences based on the current lengths and the tracked numbers
    // these index the current vector which gives us the id:
    //   current[nexts[i].s0 ... nexts[i].s1)
    // are all the potential nexts for sequence i
    vector<long2_> nexts(_total, _context);

    find_nexts(current_lengths, nexts);

    // STEP 3: Tally up all our new potential sequences so we can enumerate them
    // now count up how many we think we have
    vector<long> scratch(_total, _context);

    // transform the subtraction between end and begin into our scratch
    transform(nexts.begin(), nexts.end(), scratch.begin(), pair_difference, _queue);

    // add them all up
    inclusive_scan(scratch.begin(), scratch.end(), scratch.begin(), _queue);
    // cout << "scratch:\n";
    // print(scratch);
    long found_count;
    copy_n(scratch.end() - 1, 1, &found_count, _queue);

    // std::wcout << "found: " << found_count << std::endl;
    if(found_count == 0)
        return found;

    // STEP 4: Enumerate all potential new sequences and sort the list of them by (prev,next)
    // now allocate memory to see if these are new or not
    
    // collect all our finds into prev and next arrays
    found.resize(found_count, _queue);
    collect_finds(nexts, scratch, current, found);

    return found;
}

long seqt::add_new_finds(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value) {
    // std::wcout << "finds: " << found.size() << endl;
    // grow our vectors
    _lengths.resize(_total + new_find_indices.size(), _queue);
    _seqs.resize(_total + new_find_indices.size(), _queue);
    _initial_seq_counts.resize(_total + new_find_indices.size(), _queue);
    _tracked.resize(_total + new_find_indices.size(), _queue);
    _counts.resize(_total + new_find_indices.size(), _queue);
    _initial_seq_counts.resize(_total + new_find_indices.size(), _queue);
    _initial_characters_read.resize(_total + new_find_indices.size(), _queue);

    fill(_initial_characters_read.begin() + _total, _initial_characters_read.begin() + _total + new_find_indices.size(), _characters_read, _queue);
    // initialize the new ones
    initialize_newly_found_sequences(new_find_indices, found, tracked_value);

    // TODO: this is slow an temporary
    std::vector<long2_> new_seqs(new_find_indices.size());
    copy(_seqs.begin() + _total, _seqs.end(), new_seqs.begin(), _queue);
    std::for_each(new_seqs.begin(), new_seqs.end(), [&](long2_ const & p) {
        size_t i = &p - &new_seqs[0];
        _seq_strings[_total + i] = _seq_strings[p.x] + _seq_strings[p.y];
    });

    // update our total
    _total += new_find_indices.size();

    return new_find_indices.size();
}

void seqt::remove_least_significant(long max_sequences) {
    long to_remove = _total - max_sequences;
    if(to_remove <= 0)
        return;

    // zip up the significance with the index into all the other arrays
    vector<long> dex(_total - 1, _context);
    vector<float> sig(_total - 1, _context);
    iota(dex.begin(), dex.end(), 1, _queue); // start at 1 to avoid the null sequence
    copy(_significance.begin() + 1, _significance.end(), sig.begin(), _queue);

    std::wcout << "all: " << endl;
    copy(sig.begin(), sig.end(), std::ostream_iterator<float, wchar_t>(std::wcout, L", "), _queue);
    std::wcout << endl;

    // TODO: doesn't sort as expected, and some nan's are seen
    // also adds some 0:0?
    sort_by_key(sig.begin(), sig.end(), dex.begin(), _queue);

    std::wcout << "sorted: " << endl;
    copy(sig.begin(), sig.end(), std::ostream_iterator<float, wchar_t>(std::wcout, L", "), _queue);
    std::wcout << endl;

    // find the nth least significant
    // nth_element_struct(nth.begin(), nth.begin() + to_remove, nth.end(), less_significance_index, _queue);

    std::wcout << "less: " << endl;
    copy(sig.begin(), sig.begin() + to_remove, std::ostream_iterator<float, wchar_t>(std::wcout, L", "), _queue);
    std::wcout << endl;

    std::wcout << "greater: " << endl;
    copy(sig.begin() + to_remove, sig.end(), std::ostream_iterator<float, wchar_t>(std::wcout, L", "), _queue);
    std::wcout << endl;

    // get the indices of the least significant
    vector<long> least_significant_index(to_remove, _context);
    copy(dex.begin(), dex.begin() + to_remove, least_significant_index.begin(), _queue);

    // flag all the least significant
    vector<long> flagged(_total, _context);
    fill(flagged.begin(), flagged.end(), 0, _queue);
    scatter_value(least_significant_index, 1, flagged);



    vector<long> output(_total, _context);

    //TODO: free up memory

    for(;;) {
        std::wcout << "flagged:" << endl;
        copy(flagged.begin(), flagged.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L", "), _queue);
        std::wcout << endl;

        depends_on_flagged(flagged, output);
        long updated_remove = count_if(output.begin(), output.end(), _1 == 1, _queue);

        std::wcout << "output:" << endl;
        copy(output.begin(), output.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L", "), _queue);
        std::wcout << endl;

        if(updated_remove < to_remove)
            throw std::logic_error("updated_remove < to_remove");

        if(updated_remove == to_remove)
            break;

        copy(output.begin(), output.end(), flagged.begin(), _queue);
        to_remove = updated_remove;
    }

    // all indices are now in index with [0...to_remove) containing the ones to be removed, and [to_remove+1..._total) containing the ones we will keep
    vector<long> keep_index = pack(flagged, _1 == 0);

    std::wcout << "keep_index:" << endl;
    copy(keep_index.begin(), keep_index.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L", "), _queue);
    std::wcout << endl;
    
    recreate_from_map(keep_index);

}  

void seqt::recreate_from_map(vector<long> & index) {
    // create a reverse_index
    auto index_size = index.size();

    vector<long> reverse_index(_total, _context);
    vector<long> iota_vec(index_size, _context);
    iota(iota_vec.begin(), iota_vec.end(), 0, _queue); // start with one to make room for our null sequence
    // any data to be removed will be marked with 0
    fill(reverse_index.begin(), reverse_index.end(), 0, _queue);
    scatter(iota_vec.begin(), iota_vec.end(), index.begin(), reverse_index.begin(), _queue);

    std::wcout << "reerse_index: " << endl;
    copy(reverse_index.begin(), reverse_index.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L","), _queue);
    std::wcout << endl;

    // now reverse_index[i] is the location of sequence `i` in our new arrays
    vector<long> new_counts(index_size, _context);
    vector<long> new_lengths(index_size, _context);
    vector<long2_> new_seqs(index_size, _context);
    vector<long2_> new_initial_seq_counts(index_size, _context);
    vector<long> new_initial_characters_read(index_size, _context);
    vector<float> new_expected_counts(index_size, _context);
    vector<float> new_stddev_counts(index_size, _context);
    vector<float> new_significance(index_size, _context);
    vector<long> new_tracked(index_size, _context);

    // get the null sequence
    // zero will always be there
    // copy_n(_counts.begin(), 1, new_counts.begin(), _queue);
    // copy_n(_lengths.begin(), 1, new_lengths.begin(), _queue);
    // copy_n(_seqs.begin(), 1, new_seqs.begin(), _queue);
    // copy_n(_initial_seq_counts.begin(), 1, new_initial_seq_counts.begin(), _queue);
    // copy_n(_initial_characters_read.begin(), 1, new_initial_characters_read.begin(), _queue);
    // copy_n(_expected_counts.begin(), 1, new_expected_counts.begin(), _queue);
    // copy_n(_stddev_counts.begin(), 1, new_stddev_counts.begin(), _queue);
    // copy_n(_significance.begin(), 1, new_significance.begin(), _queue);
    // copy_n(_tracked.begin(), 1, new_tracked.begin(), _queue);

    // then gather the rest
    std::wcout << "counts: " << endl;
    copy(_counts.begin(), _counts.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L","), _queue);
    std::wcout << endl;

    std::wcout << "index: " << endl;
    copy(index.begin(), index.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L","), _queue);
    std::wcout << endl;


    gather(index.begin(), index.end(), _counts.begin(), new_counts.begin(), _queue);


    std::wcout << "new_counts: " << endl;
    copy(new_counts.begin(), new_counts.end(), std::ostream_iterator<long, wchar_t>(std::wcout, L","), _queue);
    std::wcout << endl;

    gather(index.begin(), index.end(), _lengths.begin(), new_lengths.begin(), _queue);
    gather(index.begin(), index.end(), _initial_seq_counts.begin(), new_initial_seq_counts.begin(), _queue);
    gather(index.begin(), index.end(), _initial_characters_read.begin(), new_initial_characters_read.begin(), _queue);
    gather(index.begin(), index.end(), _expected_counts.begin(), new_expected_counts.begin(), _queue);
    gather(index.begin(), index.end(), _stddev_counts.begin(), new_stddev_counts.begin(), _queue);
    gather(index.begin(), index.end(), _significance.begin(), new_significance.begin(), _queue);
    gather(index.begin(), index.end(), _tracked.begin(), new_tracked.begin(), _queue);

    // carefully rewrite the seqs
    gather_seqs(index, reverse_index, new_seqs);

    // now overwrite our internal accounting vectors
    _counts = new_counts;
    _lengths = new_lengths;
    _seqs = new_seqs;
    _initial_seq_counts = new_initial_seq_counts;
    _initial_characters_read = new_initial_characters_read;
    _expected_counts = new_expected_counts;
    _stddev_counts = new_stddev_counts;
    _significance = new_significance;
    _tracked = new_tracked;

    // and finally shrink our total
    _total = index_size;

    // now go through our indices and rewrite the atoms
    map<wchar_t, long> new_char_index;
    map<long, wchar_t> new_index_char;
    map<long, std::wstring> new_seq_strings;

    // TODO: paralelize this
    // TODO: there is something wrong here. The strings aren't being selected and saved
    // 
    std::vector<long> cpu_reverse_index(reverse_index.size());
    copy(reverse_index.begin(), reverse_index.end(), cpu_reverse_index.begin(), _queue);

    std::for_each(_char_index.begin(), _char_index.end(), [&](auto const & p) {
        long rev = cpu_reverse_index[p.second];
        if(rev == 0 && p.second != 0)
            return; // forgetting this one!
        
        new_char_index[p.first] = rev;
        new_index_char[rev] = p.first;
    });

    std::for_each(_seq_strings.begin(), _seq_strings.end(), [&](auto const & p) {
        long rev = cpu_reverse_index[p.first];
        if(rev == 0 && p.first != 0)
            return;

        new_seq_strings[rev] = _seq_strings[p.first];
    });

    _char_index = new_char_index;
    _index_char = new_index_char;
    _seq_strings = new_seq_strings;
}


void seqt::scatter_value(vector<long> & indices, long value, vector<long> & output) {
    long local_size = _scatter_value_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(indices.size(), local_size);

    _scatter_value_kernel.set_arg(0, indices);
    _scatter_value_kernel.set_arg(1, (long)indices.size());
    _scatter_value_kernel.set_arg(2, value);
    _scatter_value_kernel.set_arg(3, output);

    event e = _queue.enqueue_1d_range_kernel(_scatter_value_kernel, 0, operational_size, local_size);
    e.wait();
}

void seqt::collect_finds(vector<long2_> & nexts, vector<long> & scratch, vector<long> & current, vector<long2_> & found) {
    long local_size = _collect_finds_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(found.size(), local_size);

    _collect_finds_kernel.set_arg(0, nexts);
    _collect_finds_kernel.set_arg(1, scratch);
    _collect_finds_kernel.set_arg(2, (long)scratch.size());
    _collect_finds_kernel.set_arg(3, current);
    _collect_finds_kernel.set_arg(4, found);
    _collect_finds_kernel.set_arg(5, (long)found.size());

    event e = _queue.enqueue_1d_range_kernel(_collect_finds_kernel, 0, operational_size, local_size);
    e.wait();

}
void seqt::mark_exists(vector<long2_> & found, vector<long> & scratch) {
    long local_size = _mark_exists_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    _mark_exists_kernel.set_arg(0, _seqs);
    _mark_exists_kernel.set_arg(1, _total);
    _mark_exists_kernel.set_arg(2, found);
    _mark_exists_kernel.set_arg(3, scratch);
    _mark_exists_kernel.set_arg(4, (long)found.size());

    event e = _queue.enqueue_1d_range_kernel(_mark_exists_kernel, 0, operational_size, local_size);
    e.wait();
}
void seqt::initialize_newly_found_sequences(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value) {
    long local_size = _initialize_newly_found_sequences_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(new_find_indices.size(), local_size);

    _initialize_newly_found_sequences_kernel.set_arg(0, new_find_indices);
    _initialize_newly_found_sequences_kernel.set_arg(1, (long)new_find_indices.size());
    _initialize_newly_found_sequences_kernel.set_arg(2, found);
    _initialize_newly_found_sequences_kernel.set_arg(3, _lengths);
    _initialize_newly_found_sequences_kernel.set_arg(4, _seqs);
    _initialize_newly_found_sequences_kernel.set_arg(5, _initial_seq_counts);
    _initialize_newly_found_sequences_kernel.set_arg(6, _tracked);
    _initialize_newly_found_sequences_kernel.set_arg(7, _counts);
    _initialize_newly_found_sequences_kernel.set_arg(8, tracked_value);
    _initialize_newly_found_sequences_kernel.set_arg(9, _total);

    event e = _queue.enqueue_1d_range_kernel(_initialize_newly_found_sequences_kernel, 0, operational_size, local_size);
    e.wait();
}


void seqt::find_nexts(vector<long> & sorted_lengths, vector<long2_> & nexts) {
    // calculate local and operational sizes
    long local_size = _find_nexts_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    // set parameters
    _find_nexts_kernel.set_arg(0, _tracked);
    _find_nexts_kernel.set_arg(1, _total);
    _find_nexts_kernel.set_arg(2, sorted_lengths);
    _find_nexts_kernel.set_arg(3, (long)sorted_lengths.size());
    _find_nexts_kernel.set_arg(4, nexts);

    // run the kernel
    event e = _queue.enqueue_1d_range_kernel(_find_nexts_kernel, 0, operational_size, local_size);
    e.wait();
}

vector<long2_> seqt::make_pair_constant_second(vector<long> & first, long second) {
    long local_size = _make_pair_constant_second_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    vector<long2_> output(first.size(), _context);

    _make_pair_constant_second_kernel.set_arg(0, first);
    _make_pair_constant_second_kernel.set_arg(1, second);
    _make_pair_constant_second_kernel.set_arg(2, output);
    _make_pair_constant_second_kernel.set_arg(3, first.size());

    event e = _queue.enqueue_1d_range_kernel(_make_pair_constant_second_kernel, 0, operational_size, local_size);
    e.wait();

    return output;
}

void seqt::is_sequence_significant(vector<long2_> & seq, float sigma, long min_count, vector<long> & output) {
    long local_size = _is_sequence_significant_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(seq.size(), local_size);

    _is_sequence_significant_kernel.set_arg(0, seq);
    _is_sequence_significant_kernel.set_arg(1, seq.size());
    _is_sequence_significant_kernel.set_arg(2, _significance);
    _is_sequence_significant_kernel.set_arg(3, _counts);
    _is_sequence_significant_kernel.set_arg(4, _total);
    _is_sequence_significant_kernel.set_arg(5, sigma);
    _is_sequence_significant_kernel.set_arg(6, min_count);
    _is_sequence_significant_kernel.set_arg(7, output);

    event e = _queue.enqueue_1d_range_kernel(_is_sequence_significant_kernel, 0, operational_size, local_size);
    e.wait();
}



void seqt::calculate_stats() {
    long local_size = _calculate_stats_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    _expected_counts = vector<float>(_total, _context);
    _stddev_counts = vector<float>(_total, _context);
    _significance = vector<float>(_total, _context);

    _calculate_stats_kernel.set_arg(0, _seqs);
    _calculate_stats_kernel.set_arg(1, _initial_seq_counts);
    _calculate_stats_kernel.set_arg(2, _counts);
    _calculate_stats_kernel.set_arg(3, _initial_characters_read);
    _calculate_stats_kernel.set_arg(4, _lengths);
    _calculate_stats_kernel.set_arg(5, _characters_read);
    _calculate_stats_kernel.set_arg(6, _stddev_counts);
    _calculate_stats_kernel.set_arg(7, _expected_counts);
    _calculate_stats_kernel.set_arg(8, _significance);
    _calculate_stats_kernel.set_arg(9, _total);

    event e = _queue.enqueue_1d_range_kernel(_calculate_stats_kernel, 0, operational_size, local_size);
    e.wait();

    //DEBUG: looking for source of nan
    vector<int> check(_total, _context);
    transform(_significance.begin(), _significance.end(), check.begin(), check_nan, _queue);

    auto found = accumulate(check.begin(), check.end(), 0, _queue);
    if(found > 0) {
        wcout << "found nan:\n";
        copy(_significance.begin(), _significance.end(), ostream_iterator<float, wchar_t>(wcout, L", "), _queue);
        wcout << endl;
        wcout << endl;
    }
}

void seqt::gather_seqs(vector<long> & index, vector<long> & reverse_index, vector<long2_> & new_seqs) {
    long local_size = _gather_seqs_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(index.size(), local_size);

    _gather_seqs_kernel.set_arg(0, index);
    _gather_seqs_kernel.set_arg(1, (long)index.size());
    _gather_seqs_kernel.set_arg(2, _seqs);
    _gather_seqs_kernel.set_arg(3, reverse_index);
    _gather_seqs_kernel.set_arg(4, new_seqs);

    event e = _queue.enqueue_1d_range_kernel(_gather_seqs_kernel, 0, operational_size, local_size);
    e.wait();
}

void seqt::depends_on_flagged(vector<long> & flagged, vector<long> & output) {
    long local_size = _depends_on_flagged.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    _depends_on_flagged.set_arg(0, _seqs);
    _depends_on_flagged.set_arg(1, _total);
    _depends_on_flagged.set_arg(2, flagged);
    _depends_on_flagged.set_arg(3, output);

    event e = _queue.enqueue_1d_range_kernel(_depends_on_flagged, 0, operational_size, local_size);
    e.wait();
}

void seqt::depends_on_sorted_list(buffer_iterator<long> sorted_begin, buffer_iterator<long> sorted_end, 
                            buffer_iterator<long> begin, buffer_iterator<long> end,
                            vector<long> & output)
{
    auto total = end - begin;
    auto sorted_total = sorted_end - sorted_begin;

    long local_size = _depends_on_sorted_list_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(total, local_size);


    // TODO: fix error in code below:
    //    main: /usr/include/boost/compute/buffer.hpp:173: buffer boost::compute::buffer::create_subbuffer(cl_mem_flags, size_t, size_t): Assertion `origin % (get_context(). get_device(). get_info<0x1019>() / 8) == 0' failed.
    // buffer sorted_buffer = sorted_begin.get_buffer();
    // buffer data_buffer = begin.get_buffer();
    // buffer output_buffer = output_begin.get_buffer();

    // auto sorted_subbuffer = sorted_buffer.create_subbuffer(CL_MEM_READ_ONLY, sorted_begin.get_index(), sorted_total);
    // auto data_subbuffer = data_buffer.create_subbuffer(CL_MEM_READ_ONLY, begin.get_index(), total);
    // auto output_subbuffer = output_buffer.create_subbuffer(CL_MEM_READ_WRITE, output_begin.get_index(), output_buffer.size() - output_begin.get_index());

    // TODO: figure out how to do this with subbuffers or iterators to avoid these copies.  The code above errors due to alignment in OpenCL
    vector<long> sorted_subbuffer(sorted_begin, sorted_end, _queue);
    vector<long> data_subbuffer(begin, end, _queue);

    _depends_on_sorted_list_kernel.set_arg(0, _seqs);
    _depends_on_sorted_list_kernel.set_arg(1, sorted_subbuffer);
    _depends_on_sorted_list_kernel.set_arg(2, (long)sorted_total);
    _depends_on_sorted_list_kernel.set_arg(3, data_subbuffer);
    _depends_on_sorted_list_kernel.set_arg(4, total);
    _depends_on_sorted_list_kernel.set_arg(5, output);

    event e = _queue.enqueue_1d_range_kernel(_depends_on_sorted_list_kernel, 0, operational_size, local_size);
    e.wait();
}

vector<long> seqt::pack(vector<long> & data, function<long(long)> pred) {
    // allocate a scratch vector
    vector<long> scratch(data.size(), _context);

    // transform the predicate into our scratch
    transform(data.begin(), data.end(), scratch.begin(), pred, _queue);

    // scan the scratch data to count and prep for packing
    inclusive_scan(scratch.begin(), scratch.end(), scratch.begin(), _queue);

    // get the count from the predicate
    long count;
    copy_n(scratch.end() - 1, 1, &count, _queue);

    // allocate space for the packed indices
    vector<long> packed(count, _context);

    // calculate local and operational sizes
    long local_size = _pack_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(scratch.size(), local_size);

    // set parameters
    _pack_kernel.set_arg(0, long_(scratch.size()));
    _pack_kernel.set_arg(1, scratch);
    _pack_kernel.set_arg(2, packed);

    // run our pack kernel
    event e = _queue.enqueue_1d_range_kernel(_pack_kernel, 0, operational_size, local_size);
    e.wait();

    return packed;
}


template<typename T>
void seqt::print(std::wostream & os, vector<T> & v) {
    std::vector<T> d(v.size());
    copy(v.begin(), v.end(), d.begin(), _queue);
    
    const wchar_t delim[] = {',', '\0'};
    std::copy(d.begin(), d.end(), ostream_iterator<T, wchar_t>(os, delim));
    os << endl;
}


long calc_operational_size(long global_size, long local_size) {
    if (global_size % local_size == 0) {
        return global_size;
    }
    return local_size * (1 + global_size / local_size);
}


seqt::seqt() :
    _device(system::default_device()),
    _context(_device),
    _queue(_context, _device),
    _program(program::build_with_source_file("cl/kernels.cl", _context)),
    _counts(0, _context),
    _lengths(0, _context),
    _seqs(0, _context),
    _initial_seq_counts(0, _context),
    _initial_characters_read(0, _context),
    _tracked(0, _context),
    _expected_counts(0, _context),
    _stddev_counts(0, _context),
    _total(0),
    _characters_read(0),
    _sigma(5.),
    _min_sigma(2.),
    _min_occurances(5),
    _max_sequences_tracked(1e3)
{
    _program.build();

    _pack_kernel = _program.create_kernel("pack");
    _scatter_value_kernel = _program.create_kernel("scatter_value");
    _find_nexts_kernel = _program.create_kernel("find_nexts");
    _collect_finds_kernel = _program.create_kernel("collect_finds");
    _mark_exists_kernel = _program.create_kernel("mark_exists");
    _initialize_newly_found_sequences_kernel = _program.create_kernel("initialize_newly_found_sequences");
    _make_pair_constant_second_kernel = _program.create_kernel("make_pair_constant_second");
    _calculate_stats_kernel = _program.create_kernel("calculate_stats");
    _is_sequence_significant_kernel = _program.create_kernel("is_sequence_significant");
    _depends_on_sorted_list_kernel = _program.create_kernel("depends_on_sorted_list");
    _gather_seqs_kernel = _program.create_kernel("gather_seqs");
    _depends_on_flagged = _program.create_kernel("depends_on_flagged");

    get_char_index(0);

    // print(std::wcout, _tracked);
    // HACK: set the 0 symbol to be tracked behind any other sequence
    fill_n(_tracked.begin(), 1, -1, _queue);
    // print(std::wcout, _tracked);
}

// TODO: this is filling up the stack because of the recursion
std::wstring seqt::write_by_index(long i, std::vector<long2_> const & seqs, 
    std::map<long, std::wstring> & cache)
{
    std::wstring str;

    // index 0 is always the empty string
    if(i == 0) {
        str.resize(1);
        str[0] = ' ';
        return str;
    }
    
    auto j = cache.find(i);
    if(j != cache.end()) {
        return j->second;
    }

    // atom
    if(seqs[i] == long2_(0,0)) {
        wchar_t c = _index_char[i];
        str.resize(1);
        str[0] = c;

        // std::wcout << "atom: " << i << " " << str << endl;
        cache[i] = str;
        return str;
    }

    // std::wstring p, n;

    // // TODO: make this stack based instead of recursive later
    // p = write_by_index(seqs[i].x, seqs, cache);
    // n = write_by_index(seqs[i].y, seqs, cache);

    // str = p + n;

    str = _seq_strings[i];

    cache[i] = str;
    // std::wcout << "found: " << i << " " << str << endl;
    return str;
}

void seqt::print_by_index(std::wostream & os, long i, std::vector<long2_> const & seqs,
    std::vector<long> const & counts, std::vector<float> & sig, std::map<long, std::wstring> & cache) 
{
    std::wstring str;

    auto j = cache.find(i);

    if(j == cache.end()) {
        str = write_by_index(i, seqs, cache);
        // cache.insert({i, str}); // shouldn't be needed....
    } else {
        str = j->second;
    }


    os << "\"" << str << "\": " << counts[i] << " " << sig[i] << "\n";
}

void seqt::print_all(std::wostream & os) {

    os << "counts: " << _counts.size() << endl;
    print(os, _counts);

    std::vector<long2_> seqs(_total);
    std::vector<long> counts(_total);
    std::vector<float> sig(_total);

    copy(_seqs.begin(), _seqs.end(), seqs.begin(), _queue);
    copy(_counts.begin(), _counts.end(), counts.begin(), _queue);
    copy(_significance.begin(), _significance.end(), sig.begin(), _queue);

    std::map<long, std::wstring> cache;

    for(long i = 0; i < _total; i++) {
        print_by_index(os, i, seqs, counts, sig, cache);
    }
    os << std::endl;
}
