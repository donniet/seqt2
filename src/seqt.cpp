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
        _nexts.resize(_total, _queue);
        _prevs.resize(_total, _queue);
        _tracked.resize(_total, _queue);

        fill_n(_lengths.begin() + index, 1, 1, _queue);
        fill_n(_nexts.begin() + index, 1, 0, _queue);
        fill_n(_prevs.begin() + index, 1, 0, _queue);

        fill_n(_counts.begin() + index, 1, 1, _queue);
        fill_n(_tracked.begin() + index, 1, 1, _queue); // initialize tracked to a >0 value
    } else {
        index = i->second;// increment the count

        long c;
        copy_n(_counts.begin() + index, 1, &c, _queue);
        fill_n(_counts.begin() + index, 1, c + 1, _queue);
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

void seqt::read(wchar_t c) {

    std::wcout << "reading: " << c << endl;

    vector<long> current(0, _context);
    vector<long> current_lengths(0, _context);
    vector<long> nexts_begin(0, _context);
    vector<long> nexts_end(0, _context);
    vector<long> scratch(0, _context);
    vector<long> new_find_indices(0, _context);
    vector<long2_> found(0, _context);
    long found_count;

    // cout << "lengths:\n";
    // print(_lengths);

    // find the index of the current character
    long index = get_char_index(c); 

    // mark all the currently tracked=0 sequences as followed by index
    // this will catch all the repeat sequences
    current = pack(_tracked, _1 == 0);

    // make a pair from this and our index
    found = make_pair_constant_second(current, index);

    
    // use a tracked_value of 1 here because we will decrement all tracked values afterward
    increment_and_add_new_finds(found, 1);

    // subtract one from each tracked reference
    transform(_tracked.begin(), _tracked.end(), _tracked.begin(), _1 - 1, _queue);

    // set the currently read char as tracked=0
    fill_n(_tracked.begin() + index, 1, 0, _queue);


    for(;;) {
        // cout << "tracked:\n";
        // print(_tracked);
        // STEP 0: get all the current sequences
        current = pack(_tracked, _1 == 0);

        // mark all the current ones to a positive number to flag them as accounted for
        transform(_tracked.begin(), _tracked.end(), _tracked.begin(), flag_tracked, _queue);

        // STEP 1: sort all current sequences by length
        // allocate space for the current lengths
        current_lengths = vector<long>(current.size(), _context);

        // gather the lengths of the current items into our allocated space
        gather(current.begin(), current.end(), _lengths.begin(), current_lengths.begin(), _queue);

        // sort the current indices by the lengths
        sort_by_key(current_lengths.begin(), current_lengths.end(), current.begin(), _queue);

        // cout << "prevs and nexts:\n";
        // print(_prevs);
        // print(_nexts);

        // cout << "current:\n";
        // print(current);
        // print(current_lengths);

        // STEP 2: Find potentially new sequences by looking how far back we saw
        //   a given sequence and marking in our sorted current lengths the possible
        //   subsequent sequences
        // now go through and find new sequences based on the current lengths and the tracked numbers
        // these index the current vector which gives us the id:
        //   current[nexts_begin[i] ... nexts_end[i])
        // are all the potential nexts for sequence i
        nexts_begin = vector<long>(_total, _context);
        nexts_end = vector<long>(_total, _context);


        find_nexts(current_lengths, nexts_begin, nexts_end);
        // cout << "nexts:" << "\n";
        // print(nexts_begin);
        // print(nexts_end);


        // STEP 3: Tally up all our new potential sequences so we can enumerate them
        // now count up how many we think we have
        scratch = vector<long>(_total, _context);

        // transform the subtraction between end and begin into our scratch
        transform(nexts_begin.begin(), nexts_begin.end(), nexts_end.begin(), scratch.begin(), _2 - _1, _queue);

        // add them all up
        inclusive_scan(scratch.begin(), scratch.end(), scratch.begin(), _queue);
        // cout << "scratch:\n";
        // print(scratch);
        copy_n(scratch.end() - 1, 1, &found_count, _queue);

        std::cout << "found: " << found_count << std::endl;
        if(found_count == 0)
            break;

        // STEP 4: Enumerate all potential new sequences and sort the list of them by (prev,next)
        // now allocate memory to see if these are new or not
        found.resize(found_count, _queue);
        

        // collect all our finds into prev and next arrays
        collect_finds(nexts_begin, scratch, current, found);

        // sort it to use it as an index
        sort(found.begin(), found.end(), long2_compare, _queue);

        // STEP 5: Lookup each of our existing sequences in this list of potential new ones
        if(increment_and_add_new_finds(found, 0) == 0)
            break;
    }

    // now set all the flaged ones back to zero before the next char is read in
    transform(_tracked.begin(), _tracked.end(), _tracked.begin(), unflag_tracked, _queue);

    // look at the sequences we are tracking and see if any of them are followed by this character
    // mark those as tracked and increment their count.

    // of all the sequences we are tracking that aren't followed by this character, create new sequences

}

std::wstring seqt::write_by_index(long i, std::vector<long> const & prev, std::vector<long> const & next,
    std::map<long, std::wstring> & cache)
{
    std::wstring str;

    // index 0 is always the empty string
    if(i == 0) 
        return str;
    
    auto j = cache.find(i);
    if(j != cache.end()) {
        return j->second;
    }

    // atom
    if(prev[i] == 0 && next[i] == 0) {
        wchar_t c = _index_char[i];
        str.resize(1);
        str[0] = c;
        cache[i] = str;
        return str;
    }

    std::wstring p, n;

    // TODO: make this stack based instead of recursive later
    p = write_by_index(prev[i], prev, next, cache);
    n = write_by_index(next[i], prev, next, cache);

    str = p + n;
    cache[i] = str;
    return str;
}

void seqt::print_by_index(long i, std::vector<long> const & prev, std::vector<long> const & next,
    std::map<long, std::wstring> & cache) 
{
    std::wstring str;

    auto j = cache.find(i);

    if(j == cache.end()) {
        std::wstring str = write_by_index(i, prev, next, cache);
        // cache.insert({i, str}); // shouldn't be needed....
    } else {
        str = j->second;
    }

    std::wcout << "\"" << str << "\"" << "\n";
}

void seqt::print_all() {
    print(_prevs);
    print(_nexts);

    std::vector<long> prev(_total);
    std::vector<long> next(_total);

    copy(_prevs.begin(), _prevs.end(), prev.begin(), _queue);
    copy(_nexts.begin(), _nexts.end(), next.begin(), _queue);

    std::map<long, std::wstring> cache;

    for(long i = 0; i < _total; i++) {
        print_by_index(i, prev, next, cache);
    }
    std::cout << std::endl;
}

long seqt::increment_and_add_new_finds(vector<long2_> & found, long tracked_value) {
    if(found.size() == 0) 
        return 0;

    vector<long> scratch(found.size(), _context);
    //   mark the ones that are new, and increment counts of the ones we alread have
    // resize our scratch vector to see if we already have any of these finds
    fill(scratch.begin(), scratch.end(), 1, _queue);

    // mark the indices in scratch for new ones, and increment the count for existing ones
    mark_new_or_increment_count(found, scratch);

    // cout << "found and scratch after mark:\n";
    // print(found);
    // print(scratch);

    // STEP 6: Collect all the new sequences and add them to our list to be tracked
    // now we have indices to all new sequences and existing ones
    vector<long> new_find_indices = pack(scratch, _1 == 1);

    if(new_find_indices.size() == 0)
        return 0;

    // grow our vectors
    _lengths.resize(_total + new_find_indices.size(), _queue);
    _nexts.resize(_total + new_find_indices.size(), _queue);
    _prevs.resize(_total + new_find_indices.size(), _queue);
    _tracked.resize(_total + new_find_indices.size(), _queue);

    // initialize the new ones
    initialize_newly_found_sequences(new_find_indices, found, tracked_value);

    // update our total
    _total += new_find_indices.size();

    return new_find_indices.size();
}

void seqt::collect_finds(vector<long> & nexts_begin, vector<long> & scratch, vector<long> & current, vector<long2_> & found) {
    long local_size = _collect_finds_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(found.size(), local_size);

    _collect_finds_kernel.set_arg(0, nexts_begin);
    _collect_finds_kernel.set_arg(1, scratch);
    _collect_finds_kernel.set_arg(2, (long)scratch.size());
    _collect_finds_kernel.set_arg(3, current);
    _collect_finds_kernel.set_arg(4, found);
    _collect_finds_kernel.set_arg(5, (long)found.size());

    event e = _queue.enqueue_1d_range_kernel(_collect_finds_kernel, 0, operational_size, local_size);
    e.wait();

}
void seqt::mark_new_or_increment_count(vector<long2_> & found, vector<long> & scratch) {
    long local_size = _mark_new_or_increment_count_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    _mark_new_or_increment_count_kernel.set_arg(0, _prevs);
    _mark_new_or_increment_count_kernel.set_arg(1, _nexts);
    _mark_new_or_increment_count_kernel.set_arg(2, _counts);
    _mark_new_or_increment_count_kernel.set_arg(3, _total);
    _mark_new_or_increment_count_kernel.set_arg(4, found);
    _mark_new_or_increment_count_kernel.set_arg(5, scratch);
    _mark_new_or_increment_count_kernel.set_arg(6, (long)found.size());

    event e = _queue.enqueue_1d_range_kernel(_mark_new_or_increment_count_kernel, 0, operational_size, local_size);
    e.wait();
}
void seqt::initialize_newly_found_sequences(vector<long> & new_find_indices, vector<long2_> & found, long tracked_value) {
    long local_size = _initialize_newly_found_sequences_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(new_find_indices.size(), local_size);

    _initialize_newly_found_sequences_kernel.set_arg(0, new_find_indices);
    _initialize_newly_found_sequences_kernel.set_arg(1, (long)new_find_indices.size());
    _initialize_newly_found_sequences_kernel.set_arg(2, found);
    _initialize_newly_found_sequences_kernel.set_arg(3, _lengths);
    _initialize_newly_found_sequences_kernel.set_arg(4, _prevs);
    _initialize_newly_found_sequences_kernel.set_arg(5, _nexts);
    _initialize_newly_found_sequences_kernel.set_arg(6, _tracked);
    _initialize_newly_found_sequences_kernel.set_arg(7, tracked_value);
    _initialize_newly_found_sequences_kernel.set_arg(8, _total);

    event e = _queue.enqueue_1d_range_kernel(_initialize_newly_found_sequences_kernel, 0, operational_size, local_size);
    e.wait();
}


void seqt::find_nexts(vector<long> & sorted_lengths, vector<long> & nexts_begin, vector<long> & nexts_end) {
    // calculate local and operational sizes
    long local_size = _find_nexts_kernel.get_work_group_info<long>(_device, CL_KERNEL_WORK_GROUP_SIZE);
    long operational_size = calc_operational_size(_total, local_size);

    // set parameters
    _find_nexts_kernel.set_arg(0, _tracked);
    _find_nexts_kernel.set_arg(1, _total);
    _find_nexts_kernel.set_arg(2, sorted_lengths);
    _find_nexts_kernel.set_arg(3, (long)sorted_lengths.size());
    _find_nexts_kernel.set_arg(4, nexts_begin);
    _find_nexts_kernel.set_arg(5, nexts_end);

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
void seqt::print(vector<T> & v) {
    std::vector<T> d(v.size());
    copy(v.begin(), v.end(), d.begin(), _queue);
    
    std::copy(d.begin(), d.end(), ostream_iterator<T>(cout, ","));
    cout << endl;
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
    _nexts(0, _context),
    _prevs(0, _context),
    _tracked(0, _context),
    _total(0)
{
    _program.build();

    _pack_kernel = _program.create_kernel("pack");
    _find_nexts_kernel = _program.create_kernel("find_nexts");
    _collect_finds_kernel = _program.create_kernel("collect_finds");
    _mark_new_or_increment_count_kernel = _program.create_kernel("mark_new_or_increment_count");
    _initialize_newly_found_sequences_kernel = _program.create_kernel("initialize_newly_found_sequences");
    _make_pair_constant_second_kernel = _program.create_kernel("make_pair_constant_second");

    get_char_index(0);
}