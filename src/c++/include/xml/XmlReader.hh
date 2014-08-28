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
 ** \file XmlReader.hh
 **
 ** Helper classes for parsing xml
 **
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_XML_XML_READER_HH
#define iSAAC_XML_XML_READER_HH

#include <libxml/xmlreader.h>

#include <iostream>
#include <boost/noncopyable.hpp>

#include "common/Debug.hh"
#include "common/Exceptions.hh"

namespace isaac
{
namespace xml
{

class XmlReaderException : public common::IsaacException
{
public:
    XmlReaderException(const std::string &message) : common::IsaacException(message){}
};

class XmlReader: boost::noncopyable
{
    std::istream &is_;
    xmlTextReaderPtr xmlReader_;
    bool good_;
    static int xmlInputReadCallback(void * context, char * buffer, int len);
    static int xmlInputCloseCallback(void * context);
public:
    XmlReader(std::istream &is);
    ~XmlReader();

    std::string getCurrentDebugContext() const;

    bool good() const {return good_;}
    void clear() {good_ = true;}

    bool read();
    const char *getName() const;
    bool assertName(const char *name, const bool noThrow = false) const;
    bool checkName(const char *name) const;
    xmlReaderTypes getNodeType() const;

    class ElementText
    {
        const char *value_;
    public:
        ElementText(const char *value) : value_(value){}
        ElementText(const xmlChar *value) : value_(reinterpret_cast<const char*>(value)){}

        std::string string() const
        {
            return value_;
        }

        operator unsigned() const
        {
            return boost::lexical_cast<unsigned>(value_);
        }

        operator unsigned long() const
        {
            return boost::lexical_cast<unsigned long>(value_);
        }

        operator double() const
        {
            return boost::lexical_cast<double>(value_);
        }

        operator const char*() const
        {
            return value_;
        }
    };

    ElementText getValue() const;
    int getCurrentDepth() const;
    XmlReader &nextElement(const bool noThrow = false);
    XmlReader &nextElement(const char *name, const bool noThrow = false);
    XmlReader &nextChildElement(const char *name, const bool noThrow = false);
    XmlReader &nextElementBelowLevel(const int minDepth);

    class Attribute
    {
        xmlChar *value_;
    public:
        Attribute(xmlChar *value) : value_(value){}
        Attribute(const Attribute &that) : value_(xmlStrdup(that.value_)){}
        ~Attribute()
        {
            xmlFree(value_);
        }

        bool operator != (const char * value)
        {
            return strcmp(reinterpret_cast<const char*>(value_), value);
        }

        Attribute & operator = (Attribute that)
        {
            std::swap(value_, that.value_);
            return *this;
        }

        template <typename T> operator T() const
        {
            return boost::lexical_cast<T>(value_);
        }
    };

    Attribute getAttribute(const char *name) const;

    template <typename DefaulT>
    DefaulT getAttribute(const char *name, const DefaulT &defaultValue) const
    {
        xmlChar *attributeValue = getAttributePointer(name);
        if (!attributeValue)
        {
            return defaultValue;
        }

        Attribute ret(attributeValue);
        return ret;
    }

    ElementText readElementText();

    /**
     * \brief throws if currentName does not match
     */
    XmlReader &operator()(const char *name)
    {
        assertName(name);
        return *this;
    }

    Attribute operator[](const char *name) const
    {
        return getAttribute(name);
    }

    /**
     * \brief advances to next element and throws if fails
     */
    XmlReader & operator++(int)
    {
        return nextElement(false);
    }

    /**
     * \brief advances to next element and does not throw if fails
     */
    XmlReader & operator++()
    {
        return nextElement(true);
    }

    XmlReader & operator+=(const char *nextElementName)
    {
        return nextElement(nextElementName, false);
    }

    operator bool() const
    {
        return good();
    }

    bool operator !() const
    {
        return !good();
    }
private:
    std::string getXmlNodePath(xmlNodePtr currentNode) const;
    xmlChar *getAttributePointer(const char *name) const;
};


} // namespace xml
} // namespace isaac

#endif // #ifndef iSAAC_XML_XML_READER_HH
