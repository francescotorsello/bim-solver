/**
 *  @file      trackUsedTime.h
 *  @brief     Tracks the elapsed lifetime of an object.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _TRACK_USED_TIME_H_INCLUDED
#define _TRACK_USED_TIME_H_INCLUDED

#include <chrono>

/** The high resolution clock time-point.
 */
typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_hr;

/** Reports the elapsed time for the lifetime of the created object.
 */
class TrackUsedTime
{
    time_hr rt_born; //!< The realtime when the object was created

public:

    TrackUsedTime ()
    {
        rt_born = std::chrono::high_resolution_clock::now ();
    }

    ~TrackUsedTime ()
    {
        slog << "*** Elapsed: ";
        auto end = std::chrono::high_resolution_clock::now ();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                       ( end - rt_born ).count ();
        if ( elapsed < 1e3 ) {
            slog << elapsed << " ms" << std::endl;
        }
        else {
            slog << std::setprecision(3) << 1e-3 * elapsed << " s" << std::endl;
        }
    }
};

#endif // _TRACK_USED_TIME_H_INCLUDED
