/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file MD5Sum.hh
 **
 ** \brief Declaration of the common MD5Sum mechanism.
 **
 ** Simple MD5 digest calculator which wraps the standard (C) RSA 
 ** implementation.  
 **
 ** \author David Kimmel 
 **/


#ifndef iSAAC_COMMON_MD5SUM_HH 
#define iSAAC_COMMON_MD5SUM_HH 

// TODO, get rid of this include so that we can
// move the header out of the include dir
// TODO, change comment style to iSAAC
#include "common/md5.hh"
#include <string>

namespace isaac
{
namespace common 
{

class MD5Sum
{
public:
    /**
     * <b>Purpose:</b> CTOR
     */
    MD5Sum();
    
    /**
     * <b>Purpose:</b> DTOR
     */
    ~MD5Sum();

    /**
     * <b>Purpose:</b> Process input data
     * @param const unsigned char*, buffer to read
     * @param const int, buffer length to read
     */
    void read(const char* buffer, const int bufferLength);

    /**
     * <b>Purpose:</b> Get the digest
     * @return Digest that has been calculated (thus far)
     * TODO, might want to just typedef as a std::array
     */
    class Digest
    {
    public:
        unsigned char data[16];
    };
    Digest getDigest() const;

    /**
     * <b>Purpose:</b> Get a hexadecimal string from a char array
     * @param unsighed char array to read
     * @param int number of chars to read
     * @return std::string with a hexadecimal representation of the input array
     */
    static std::string toHexString(const unsigned char* buffer, const unsigned int bufferLength);

    /**
     * <b>Purpose:</b> Clear the internal state
     */
    void clear();

private:

    // Internal data
    MD5 md5_;

};

} // namespace common 
} // namespace isaac 

#endif // $ifndef iSAAC_COMMON_MD5SUM_HH 
