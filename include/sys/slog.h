/**
 *  @file      slog.h
 *  @brief     Writing both to cerr and cout simultaneously.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _SLOG_H_INCLUDED
#define _SLOG_H_INCLUDED

#include <sys/stat.h>
// #include <sys/types.h>

/** SLog enables writing both to cerr and cout simultaneously.
 */
struct SLog
{
    bool disabled;
    bool tty;

    SLog () {
        disabled = false;
        struct stat buf;
        tty = fstat( fileno( stderr ), &buf ) >= 0
              && (buf.st_mode & S_IFMT) == S_IFCHR;
    }

    SLog& operator<< ( std::ostream& (*pfun)(std::ostream&) )
    {
        if ( disabled ) return *this;
        pfun( std::cerr );
        if ( ! tty ) pfun( std::cout );
        return *this;
    }
};

template <class T>
SLog& operator<< ( SLog& st, T val )
{
    if ( st.disabled ) return st;
    std::cerr << val;
    if ( ! st.tty ) std::cout << val;
    return st;
}

static SLog slog;

#endif // _SLOG_H_INCLUDED
