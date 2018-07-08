/*
 * Minimal graph library
 * Utility functions.
 *
 * Copyright (c) 2017-2018 Sebastian Lackner
 */

#include <stdarg.h>
#include "graph.h"
#include "internal.h"

size_t resort_r(void *base, size_t nmemb, size_t size, size_t index,
                int (*compar)(const void *, const void *, void *), void *userdata)
{
    char *ptr = (char *)base + size * index;
    /* Please note that the search could also be done in logarithmic time -
     * in practice it doesn't matter much though because O(n) elements have
     * to be moved in the worst case. */
    while (index > 0 && compar(ptr - size, ptr, userdata) > 0)
    {
        SWAP(ptr - size, ptr, size);
        ptr -= size;
        index--;
    }
    while (index < nmemb - 1 && compar(ptr, ptr + size, userdata) > 0)
    {
        SWAP(ptr, ptr + size, size);
        ptr += size;
        index++;
    }
    return index;
}

void progress(const char *format, ...)
{
    static int last_len;
    int i, new_len;
    va_list args;

    va_start(args, format);
    new_len = vfprintf(stderr, format, args);
    for (i = new_len; i < last_len; i++) fputc(' ', stderr);
    fputc('\r', stderr);
    va_end(args);

    last_len = new_len;
}
