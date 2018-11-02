/**
 *  @file      paramsHolder.h
 *  @brief     Access to the configuration files.
 *  @authors   Mikica Kocic, Wilkinson & Reinsch, Press et al
 *  @copyright GNU General Public License (GPLv3).
 */

#include <map>

/** Contains `key = value` pairs that are obtained from a configuration file.
 */
class Parameters
    : std::map<std::string, std::string>
{
public:

    /** Checks whether the given key already exists.
     */
    bool exists( const std::string& key ) const {
        return count( key ) > 0;
    }

    /** Loads the configuration parameters from file.
     *  The expected format is `key = value` (where value can be empty).
     */
    Parameters( const std::string& fileName )
    {
        FILE* inf = fileName == "stdin" ? stdin
                  : fopen( fileName.c_str(), "r" );
        if( ! inf ) {
            std::cerr << "Error: Could not open " << fileName << std::endl;
            return;
        }

        char line[ 1024 ] = {0};
        while( fgets( line, sizeof(line)-1, inf ) )
        {
            char key[ 64 ] = {0};
            int eqpos = -1, value = -1;
            if ( 1 != sscanf( line, " %s %n= %n", key, &eqpos, &value ) ) {
                // empty line
            }
            else if ( eqpos < 0 || line[eqpos] != '=' ) {
                std::cerr << "Error: Not 'key = value' pair:" << std::endl << line;
            }
            else {
                char* chp = line + value;
                while( *chp ) ++chp;
                while( isspace( chp[-1] ) ) { *--chp = 0; }
                if ( exists( key ) ) {
                    std::cerr << "Warning: Redefining " << key << std::endl;
                }
                operator[]( key ) = std::string( line + value );
            }
        };

        fclose( inf );
    }

    void get( const std::string& key, Real& var, const Real defaultValue )
    {
        var = !exists( key ) ? defaultValue
            : strToReal( operator[]( key ).c_str(), NULL );
    }

    void get( const std::string& key, Int& var, const Int defaultValue )
    {
        var = !exists( key ) ? defaultValue
            : std::atoi( operator[]( key ).c_str() );
    }

    void get( const std::string& key, std::string& var, const char* defaultValue )
    {
        var = !exists( key ) ? defaultValue : operator[]( key );
    }

    void get( const std::string& key, bool& var, const bool defaultValue )
    {
        var = !exists( key ) ? defaultValue
            : std::atoi( operator[]( key ).c_str() ) != 0;
    }

    std::string get( const std::string& key, Int& var, const Int defaultValue,
              std::map<std::string, int>& enumMap )
    {
        if( !exists( key ) || enumMap.count( operator[]( key ) ) == 0 ) {
            var = defaultValue;
            return "(unknown)";
        }
        std::string verboseName = operator[]( key );
        var = enumMap[ verboseName ];
        return verboseName;
    }

    /** Dumps all the keys that were loaded.
     */
    void dump ()
    {
        for ( Parameters::iterator i = begin(); i != end(); ++i ) {
            std::cout << "(" << i->first << ")=(" << i->second << ")" << std::endl;
        }
    }
};
