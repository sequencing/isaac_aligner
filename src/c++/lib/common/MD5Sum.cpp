/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file MD5Sum.cpp
 **
 ** \brief Implementation of simple (RSA) MD5 wrapper
 **
 ** \author David Kimmel 
 **/


#include "common/MD5Sum.hh"
#include <sstream>

namespace isaac
{
namespace common
{

MD5Sum::MD5Sum()
{
    // Initialize data to init values
    clear();
}

MD5Sum::~MD5Sum()
{
    // May want to clear here so that things are "protected"
}

void MD5Sum::update(const char* buffer, const int bufferLength)
{
    // Process data buffer
    const unsigned char* bTmp = reinterpret_cast<const unsigned char*>(buffer);
    MD5Digest(&md5_, bTmp, bufferLength);   
}

MD5Sum::Digest MD5Sum::getDigest() const
{
    // Return
    Digest r;
    // Copy data to temp for final processing
    MD5 tmp = md5_;
    // Finish calculation  
    MD5Close(&tmp, r.data);
    // Return result
    return r;
}

std::string MD5Sum::toHexString(const unsigned char* buffer, const unsigned int bufferLength)
{
    // TODO, would be nice to have that BOOST::construe_cast<>() instead
    std::stringstream ss;
    for(unsigned int i = 0; i < bufferLength; ++i)
    {
        const int val = static_cast<int>(buffer[i]);
        ss << ((val < 16)?("0"):("")) << std::hex << val;
    }

    return ss.str();
}


void MD5Sum::clear()
{
    MD5Open(&md5_);
}

} // namespace common
} // namespace isaac
