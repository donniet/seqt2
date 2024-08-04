
// index of first element in sorted that is >= target
long lower_bound(long target, global long * sorted, long length) {
    long low = 0;
    long high = length - 1;

    while(low <= high) {
        long i = (low + high) / 2;

        if(sorted[i] >= target) 
            high = i - 1;
        else
            low = i + 1;
    } 

    return low;
}

// index of first element in sorted that is >= target

bool greater_equal(long2 x, long2 y) {
    if(x.s0 > y.s0)
        return true;
    if(x.s0 < y.s0)
        return false;
    if(x.s1 >= y.s1)
        return true;

    return false;
}

bool equal(long2 x, long2 y) {
    return x.s0 == y.s0 && x.s1 == y.s1;
}

long lower_bound2(long2 target, global long2 * sorted, long length) {
    long low = 0;
    long high = length - 1;

    while(low <= high) {
        long i = (low + high) / 2;

        if(greater_equal(sorted[i], target))
            high = i - 1;
        else
            low = i + 1;
    } 

    return low;
}

// index of first element in sorted that is > target
long upper_bound(long target, global long * sorted, long length) {
    long low = 0;
    long high = length;

    while(low < high) {
        long i = (low + high) / 2;

        if(target >= sorted[i]) 
            low = i + 1;
        else
            high = i;
    } 

    return low;
}

kernel void find_nexts(
    global long * tracked, 
    long total,
    global long * sorted_current_lengths,
    long current_total,
    global long * nexts_begin,
    global long * nexts_end) 
{
    const long gid = get_global_id(0);
    if(gid > total)
        return;

    // is the previous one behind by the length of the current ones?
    long search_length = -tracked[gid];

    // do a binary search through the sorted currents for search_length
    nexts_begin[gid] = lower_bound(search_length, sorted_current_lengths, current_total);
    nexts_end[gid] = upper_bound(search_length, sorted_current_lengths, current_total);
    // nexts_begin[gid] = 0;
    // nexts_end[gid] = 1;
}

kernel void collect_finds(
    global long * scratch,
    long scratch_size,
    global long * current,
    global long2 * found,
    long found_count
) {
    const long gid = get_global_id(0);
    if(gid > found_count)
        return;

    // search through the sorted scratch vector to get the index of prev

    // the upper_bound of gid will give us the subsequent 
    // because gid ranges continuously from [0,found_count)
    // and scratch is an inclusive_scan with each step being
    // the length of the next finds
    long found_end = upper_bound(gid, scratch, scratch_size);

    // now we can get the first index of this prev by going back one
    long gid_first = scratch[found_end-1];

    // found_prev is going to be the lower bound of this gid_first
    found[gid].s0 = lower_bound(gid_first, scratch, scratch_size);

    // now we can get our found_next by subtracting gid_first from gid 
    // and using it to index our current array
    found[gid].s1 = current[gid - gid_first];
}

kernel void mark_new_or_increment_count(
    global long * prevs,
    global long * nexts,
    global long * count,
    long total,
    global long2 * found,
    global long * scratch,
    long found_count
) {
    const long gid = get_global_id(0);
    if(gid > total)
        return;

    // we are looking at all our existing nodes and marking
    // scratch if this one exists in found
    long2 value;
    value.s0 = prevs[gid];
    value.s1 = nexts[gid];

    long low = lower_bound2(value, found, found_count);
    if(equal(found[low], value))
        count[gid]++;
    else
        scratch[low] = 1;
}

kernel void initialize_newly_found_sequences(
    global long * new_indices,
    long new_total,
    global long2 * found,
    global long * lengths,
    global long * prevs,
    global long * nexts,
    global long * tracked,
    long previous_total
) {
    const long gid = get_global_id(0);
    if(gid > new_total) 
        return;

    long p = found[new_indices[gid]].s0;
    long n = found[new_indices[gid]].s1;

    long new_index = previous_total + gid;
    lengths[new_index] = lengths[p] + lengths[n];
    prevs[new_index] = p;
    nexts[new_index] = n;
    tracked[new_index] = 0;
}

kernel void pack(long total, global long *scanned, global long *output) {
    const long gid = get_global_id(0);
    if (gid >= total) {
        return;
    }

    if (gid == 0) {
        if (scanned[0] != 0) {
            output[0] = 0;
        }
    } else if (scanned[gid] > scanned[gid - 1]) {
        output[scanned[gid - 1]] = gid;
    }
}
