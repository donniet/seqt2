

/* Boilerplate feature-test macros: */
#if _WIN32 || _WIN64
#define _WIN32_WINNT 0x0A00      // _WIN32_WINNT_WIN10
#define NTDDI_VERSION 0x0A000002 // NTDDI_WIN10_RS1
#include <sdkddkver.h>
#else
#define _XOPEN_SOURCE 700
#define _POSIX_C_SOURCE 200809L
#endif

#include <iostream>
#include <locale>
#include <locale.h>
#include <stdlib.h>
#include <string>

#ifndef MS_STDLIB_BUGS // Allow overriding the autodetection.
/* The Microsoft C and C++ runtime libraries that ship with Visual Studio, as
 * of 2017, have a bug that neither stdio, iostreams or wide iostreams can
 * handle Unicode input or output.  Windows needs some non-standard magic to
 * work around that.  This includes programs compiled with MinGW and Clang
 * for the win32 and win64 targets.
 *
 * NOTE TO USERS OF TDM-GCC: This code is known to break on tdm-gcc 4.9.2. As
 * a workaround, "-D MS_STDLIB_BUGS=0" will at least get it to compile, but
 * Unicode output will still not work.
 */
#if (_MSC_VER || __MINGW32__ || __MSVCRT__)
/* This code is being compiled either on MS Visual C++, or MinGW, or
 * clang++ in compatibility mode for either, or is being linked to the
 * msvcrt (Microsoft Visual C RunTime) library.
 */
#define MS_STDLIB_BUGS 1
#else
#define MS_STDLIB_BUGS 0
#endif
#endif

#if MS_STDLIB_BUGS
#include <io.h>
#include <fcntl.h>
#endif

using std::endl;
using std::istream;
using std::wcin;
using std::wcout;

void init_locale(void)
// Does magic so that wcout can work.
{
#if MS_STDLIB_BUGS
    // Windows needs a little non-standard magic.
    constexpr char cp_utf16le[] = ".1200";
    setlocale(LC_ALL, cp_utf16le);
    _setmode(_fileno(stdout), _O_WTEXT);
    _setmode(_fileno(stdin), _O_WTEXT);
#else
    // The correct locale name may vary by OS, e.g., "en_US.utf8".
    constexpr char locale_name[] = "";
    setlocale(LC_ALL, locale_name);
    std::locale::global(std::locale(locale_name));
    wcout.imbue(std::locale());
    wcin.imbue(std::locale());
#endif
}

#include "seqt.hpp"

#include <string>

int main(void)
{
    init_locale();

    wchar_t c;
    wchar_t eof = std::char_traits<wchar_t>::eof();

    seqt s;

    std::string test = "blahblahblahblah";

    for(auto c : test) {
        s.read(c);
    }

    for(;;) {
        c = wcin.get();

        if(c == eof) 
            break;
    }

    return EXIT_SUCCESS;
}