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
 ** \file XmlWriter.hh
 **
 ** Helper classes for composing xml
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_XML_XML_WRITER_HH
#define iSAAC_XML_XML_WRITER_HH

#include <libxml/xmlwriter.h>

#include <iostream>
#include <boost/noncopyable.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace xml
{

class XmlWriterException : public common::IsaacException
{
public:
    XmlWriterException(const std::string &message) : common::IsaacException(message)
    {

    }
};

class XmlWriter: boost::noncopyable
{
    std::ostream &os_;
    xmlTextWriterPtr xmlWriter_;
    static int xmlOutputWriteCallback(void * context, const char * buffer, int len);
public:
    explicit XmlWriter(std::ostream &os);
    ~XmlWriter();

    void close();
    XmlWriter &startElement(const char *name);
    XmlWriter &startElement(const std::string &name) {return startElement(name.c_str());}
    XmlWriter &endElement();
    XmlWriter &writeText(const char *text);

    template <typename T> XmlWriter &writeElement(const char *name, const T& value)
    {
        return (startElement(name) << value).endElement();
    }

    template <typename T> XmlWriter &writeElement(const std::string &name, const T& value)
    {
        return writeElement(name.c_str(), value);
    }

    template <typename T> XmlWriter &writeAttribute(const char *name, const T& value)
    {
        const std::string strValue = boost::lexical_cast<std::string>(value);
        const int written = xmlTextWriterWriteAttribute(xmlWriter_, BAD_CAST name, BAD_CAST strValue.c_str());
        if (-1 == written)
        {
            BOOST_THROW_EXCEPTION(XmlWriterException(
                std::string("xmlTextWriterWriteAttribute returned -1 for attribute name: ") + name + " value: " + strValue));
        }
        return *this;
    }

    template <typename T> XmlWriter &operator <<(const T&t)
    {
        return writeText(boost::lexical_cast<std::string>(t).c_str());
    }

    struct BlockMacroSupport
    {
        bool set_;
        BlockMacroSupport() : set_(true){}
        void reset () {set_ = false;};
        operator bool() const {return set_;}
    };
    // This is to allow macro ISAAC_XML_WRITER_ELEMENT_BLOCK to work
    operator BlockMacroSupport() const {return BlockMacroSupport();}
};

/**
 * \brief Macro for controling xml element scope. Automatically closes the element at the end of the block
 */
#define ISAAC_XML_WRITER_ELEMENT_BLOCK(writer, name) for(xml::XmlWriter::BlockMacroSupport iSaacElementBlockVariable = writer.startElement(name); iSaacElementBlockVariable; iSaacElementBlockVariable.reset(), writer.endElement())


} // namespace xml
} // namespace isaac

#endif // #ifndef iSAAC_XML_XML_WRITER_HH
