#include "seqt.hpp"


long seqt::get_char_index(wchar_t c) {
    long index;

    auto i = _char_index.find(c);
    if(i == _char_index.end()) {
        index = _total;
        _char_index.insert({c, index});
        _index_char.insert({index, c});

        _total++;

        _counts.resize(_total, _queue);
        _lengths.resize(_total, _queue);
        _seqs.resize(_total, _queue);
        _tracked.resize(_total, _queue);

        fill_n(_lengths.begin() + index, 1, 1, _queue);

        // atoms have prev[i] and next[i] == 0
        fill_n(_seqs.begin() + index, 1, long2_(0,0), _queue);

        fill_n(_counts.begin() + index, 1, 0, _queue); // initialize count to zero because we haven't tracked this yet
        fill_n(_tracked.begin() + index, 1, 1, _queue); // initialize tracked to a >0 value
    } else {
        index = i->second;// increment the count

        // don't update the count here, we'll do it in the main read loop
        // long c;
        // copy_n(_counts.begin() + index, 1, &c, _queue);
        // fill_n(_counts.begin() + index, 1, c + 1, _queue);
    }

    return index;
} 

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

BOOST_COMPUTE_FUNCTION(bool, long2_compare, (long2_ a, long2_ b), {
    if(a.x < b.x) return true;
    if(b.x < a.x) return false;
    if(a.y < b.y) return true;
    return false;
});

BOOST_COMPUTE_FUNCTION(long, pair_difference, (long2_ a), {
    return a.y - a.x;
});

void seqt::read(wchar_t c) {

    std::wcout << "reading: " << c << endl;

    vector<long> current(0, _context);
    vector<long> scratch(0, _context);
    vector<long> existing_indices(0, _context);
    vector<long> new_find_indices(0, _context);
    vector<long2_> found(0, _context);
    long found_count;
    long cnt;


    // cout << "lengths:\n";
    // print(_lengths);


    // find the index of the current character
    // print(std::wcout, _tracked);
    long index = get_char_index(c); 
    transform(_tracked.begin(), _tracked.end(), _tracked.begin(), _1 - 1, _queue);

    // set the currently read char as tracked=0
    fill_n(_tracked.begin() + index, 1, 0, _queue);
    // copy_n(_counts.begin() + index, 1, &cnt, _queue);
    // fill_n(_counts.begin() + index, 1, cnt+1, _queue);

    // flag our current sequences 
    vector<long> current_flag(_total, _context);
    transform(_tracked.begin(), _tracked.end(), current_flag.begin(), _1 == 0, _queue);


    for(;;) {
        // STEP 0: get all the current sequences
        current = pack(current_flag, _1 == 1);

        // create a list of pairs of potential sequences based on the length
        found = find_nexts_by_length(current);
        if(found.size() == 0)
            break;

        // sort it to use it as an index
        sort(found.begin(), found.end(), long2_compare, _queue);

        scratch = vector<long>(found.size(), _context);
        //   mark the ones that are new, and increment counts of the ones we alread have
        // resize our scratch vector to see if we already have any of these finds
        fill(scratch.begin(), scratch.end(), 0, _queue);

        // mark the indices in scratch for new ones, and increment the count for existing ones
        mark_exists(found, scratch);

        // STEP 6: Collect all the new sequences and add them to our list to be tracked
        // now we have indices to all new sequences and existing ones
        vector<long> new_find_indices = pack(scratch, _1 == 0);

        // std::wcout << "new found indices: " << new_find_indices.size() << endl;

        // flag any newly found sequences as current
        // how many existing indices are not already flagged as current?
        // first remove all the zeros
        auto existing_end = remove(scratch.begin(), scratch.end(), 0, _queue);
        existing_indices = vector<long>(existing_end - scratch.begin(), _context);
        copy(scratch.begin(), existing_end, existing_indices.begin(), _queue);

        // re-use scratch to see if these existing sequences have been flagged
        gather(existing_indices.begin(), existing_indices.end(), current_flag.begin(), scratch.begin(), _queue);
        long existing_flagged;
        reduce(scratch.begin(), existing_end, &existing_flagged, _queue);

        if(existing_flagged < existing_indices.size()) {
            // TODO:
            scatter_value(existing_indices, 1, current_flag);
        }


        if(new_find_indices.size() > 0) {
            add_new_finds(new_find_indices, found, 0);
            
            size_t original_size = current_flag.size();
            current_flag.resize(_total, _queue);

            // all new sequences are flaged as current
            fill(current_flag.begin() + original_size, current_flag.end(), 1, _queue); 
        } else if(existing_flagged == existing_indices.size()) {
            // we have new new indices and no existing ones flagged
            break;
        }
    }

    // increment the counts of all the flagged sequences
    std::wcout << "current: ";
    print(std::wcout, current);
    scratch = vector<long>(current.size(), _context);
    
    // increment counts for all the current sequences
    gather(current.begin(), current.end(), _counts.begin(), scratch.begin(), _queue);
    transform(scratch.begin(), scratch.end(), scratch.begin(), _1 + 1, _queue);
    scatter(scratch.begin(), scratch.end(), current.begin(), _counts.begin(), _queue);
    std::wcout << "counts: ";
    print(std::wcout, _counts);

    // mark all the current sequences as tracked=0
    scatter_value(current, 0, _tracked);
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
    //   current[nexts_begin[i] ... nexts_end[i])
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
    _tracked.resize(_total + new_find_indices.size(), _queue);
    _counts.resize(_total + new_find_indices.size(), _queue);

    // initialize the new ones
    initialize_newly_found_sequences(new_find_indices, found, tracked_value);

    // update our total
    _total += new_find_indices.size();

    return new_find_indices.size();
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
    _initialize_newly_found_sequences_kernel.set_arg(5, _tracked);
    _initialize_newly_found_sequences_kernel.set_arg(6, _counts);
    _initialize_newly_found_sequences_kernel.set_arg(7, tracked_value);
    _initialize_newly_found_sequences_kernel.set_arg(8, _total);

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
    _tracked(0, _context),
    _total(0)
{
    _program.build();

    _pack_kernel = _program.create_kernel("pack");
    _scatter_value_kernel = _program.create_kernel("scatter_value");
    _find_nexts_kernel = _program.create_kernel("find_nexts");
    _collect_finds_kernel = _program.create_kernel("collect_finds");
    _mark_exists_kernel = _program.create_kernel("mark_exists");
    _initialize_newly_found_sequences_kernel = _program.create_kernel("initialize_newly_found_sequences");
    _make_pair_constant_second_kernel = _program.create_kernel("make_pair_constant_second");

    get_char_index(0);

    print(std::wcout, _tracked);
    fill_n(_tracked.begin(), 1, 0, _queue);
    print(std::wcout, _tracked);
}


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

    std::wstring p, n;

    // TODO: make this stack based instead of recursive later
    p = write_by_index(seqs[i].x, seqs, cache);
    n = write_by_index(seqs[i].y, seqs, cache);

    str = p + n;
    cache[i] = str;
    // std::wcout << "found: " << i << " " << str << endl;
    return str;
}

void seqt::print_by_index(std::wostream & os, long i, std::vector<long2_> const & seqs,
    std::vector<long> const & counts, std::map<long, std::wstring> & cache) 
{
    std::wstring str;

    auto j = cache.find(i);

    if(j == cache.end()) {
        str = write_by_index(i, seqs, cache);
        // cache.insert({i, str}); // shouldn't be needed....
    } else {
        str = j->second;
    }


    os << "\"" << str << "\": " << counts[i] << "\n";
}

void seqt::print_all(std::wostream & os) {

    os << "counts: " << _counts.size() << endl;
    print(os, _counts);

    std::vector<long2_> seqs(_total);
    std::vector<long> counts(_total);

    copy(_seqs.begin(), _seqs.end(), seqs.begin(), _queue);
    copy(_counts.begin(), _counts.end(), counts.begin(), _queue);

    std::map<long, std::wstring> cache;

    for(long i = 0; i < _total; i++) {
        print_by_index(os, i, seqs, counts, cache);
    }
    os << std::endl;
}
