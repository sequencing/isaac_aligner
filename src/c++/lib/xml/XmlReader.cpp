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
 ** \file XmlReader.cpp
 **
 ** Helper class for working with xml
 **
 ** \author Roman Petrovski
 **/

#include "xml/XmlReader.hh"


namespace isaac
{
namespace xml
{

XmlReader::XmlReader(std::istream &is) :
    is_(is),
    xmlReader_(xmlReaderForIO(xmlInputReadCallback, xmlInputCloseCallback, &is_,
                              "", 0, XML_PARSE_DTDATTR|XML_PARSE_NOENT)),
    good_(true)
{
    if (!xmlReader_)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlReaderForIO failed"));
    }
}

XmlReader::~XmlReader()
{
    xmlFreeTextReader(xmlReader_);
}

int XmlReader::xmlInputReadCallback(void *context, char *buffer, int len)
{
    std::istream* pStream = reinterpret_cast<std::istream*>(context);
    if (!*pStream && !pStream->eof())
    {
//        ISAAC_THREAD_CERR << "stream is not open" << std::endl;
        return -1;
    }

    if (!pStream->read(buffer, len) && !pStream->eof())
    {
//        ISAAC_THREAD_CERR << "stream read failed" << std::endl;
        return -1;
    }
//    ISAAC_THREAD_CERR << "stream read " << pStream->gcount() << std::endl;
    return pStream->gcount();
}

int XmlReader::xmlInputCloseCallback(void * context)
{
    return 0;
}

XmlReader::ElementText XmlReader::readElementText()
{
    if (!good())
    {
        // don't do anything unless we're good. This forces user to analyse and clear error states
        BOOST_THROW_EXCEPTION(xml::XmlReaderException("readElementText requested on bad reader"));
    }

    int isEmpty = xmlTextReaderIsEmptyElement(xmlReader_);
    if (-1 == isEmpty)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderIsEmptyElement returned -1"));
    }

    if (isEmpty)
    {
        // return empty text for element that is empty lie this <blah/>
        return ElementText("");
    }
    else
    {
        if (!read())
        {
            good_ = false;
            BOOST_THROW_EXCEPTION(XmlReaderException("Unexpected end of the stream while looking for element text. " + getCurrentDebugContext()));
        }

        xmlReaderTypes nodeType = getNodeType();
        if (XML_READER_TYPE_END_ELEMENT == nodeType)
        {
            // return empty text for element that ends here
            return ElementText("");
        }
        else if (XML_READER_TYPE_TEXT != nodeType)
        {
            good_ = false;
            BOOST_THROW_EXCEPTION(XmlReaderException("Text requested for element which does not have text. " + getCurrentDebugContext()));
        }

        return getValue();
    }
}

XmlReader::Attribute XmlReader::getAttribute(const char *name) const
{
    xmlChar *attributeValue = getAttributePointer(name);
    if (!attributeValue)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException(std::string("Attribute ") + name + " not found" + getCurrentDebugContext()));
    }
//    ISAAC_THREAD_CERR << "attribute:" << name << " value:" << attributeValue << std::endl;

    Attribute ret(attributeValue);
    return ret;
}

xmlChar *XmlReader::getAttributePointer(const char *name) const
{
    xmlReaderTypes nodeType = getNodeType();

    if (XML_READER_TYPE_ELEMENT != nodeType)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException(std::string("Request for attribute ") + name +
                                                 " while current node type is not an element: " +
                                                 boost::lexical_cast<std::string>(nodeType) +
                                                 getCurrentDebugContext()));
    }

    return xmlTextReaderGetAttribute(xmlReader_, BAD_CAST name);
}

XmlReader &XmlReader::nextElementBelowLevel(const int minDepth)
{
    if (!good())
    {
        // don't do anything unless we're good. This forces user to analyse and clear error states
        BOOST_THROW_EXCEPTION(xml::XmlReaderException("nextElementBelowLevel requested on bad reader for minDepth:" +
                                                      boost::lexical_cast<std::string>(minDepth)));
    }
    while (read())
    {
        xmlReaderTypes nodeType = getNodeType();

        if (XML_READER_TYPE_ELEMENT == nodeType)
        {
            return *this;
        }
        else if (XML_READER_TYPE_END_ELEMENT == nodeType)
        {
//            ISAAC_THREAD_CERR << "current depth:" << getCurrentDepth() << " min depth:" << minDepth << std::endl;
            if (minDepth >= getCurrentDepth())
            {
                good_ = false;
                return *this;
            }
        }
    }
    good_ = false;
    return *this;
}

XmlReader &XmlReader::nextElement(const bool noThrow/* = false*/)
{
    return nextElement(0, noThrow);
}

/**
 * \brief advances to the next element, optionaly validates the name of the found element
 *
 * \param name if not 0, name of the next found element is asserted
 */
XmlReader &XmlReader::nextElement(const char *name, const bool noThrow/* = false*/)
{
    if (!good())
    {
        // don't do anything unless we're good. This forces user to analyse and clear error states
        BOOST_THROW_EXCEPTION(xml::XmlReaderException("nextElement requested on bad reader"));
    }
    const std::string debugContext = getCurrentDebugContext();
    while (read())
    {
        xmlReaderTypes nodeType = getNodeType();

        if (XML_READER_TYPE_ELEMENT == nodeType)
        {
            if (name)
            {
                good_ = assertName(name, noThrow);
            }
            return *this;
        }
    }
    good_ = false;
    if (!noThrow)
    {
        BOOST_THROW_EXCEPTION(xml::XmlReaderException(
            std::string("Unexpected end of xml stream while looking for next element ") + (name ? name : "") + " at " +
            debugContext));
    }
    return *this;
}

XmlReader &XmlReader::nextChildElement(const char *name, const bool noThrow/* = false*/)
{
    const int parentDepth = getCurrentDepth();
    nextElement(name, noThrow);
    if (getCurrentDepth() != parentDepth + 1)
    {
        good_ = false;
        if (!noThrow)
        {
            BOOST_THROW_EXCEPTION(xml::XmlReaderException(std::string("Child element ") + name +
                " not found for parentDepth " + boost::lexical_cast<std::string>(parentDepth) + getCurrentDebugContext()));
        }
    }
    return *this;
}


XmlReader::ElementText XmlReader::getValue() const
{
    const xmlChar *ret = xmlTextReaderConstValue(xmlReader_);
    if (!ret)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderConstValue returned 0" + getCurrentDebugContext()));
    }
    return ElementText(ret);
}

int XmlReader::getCurrentDepth() const
{
    const int ret = xmlTextReaderDepth(xmlReader_);
    if (-1 == ret)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderDepth returned -1"));
    }
    return ret;
}

std::string XmlReader::getXmlNodePath(xmlNodePtr currentNode) const
{
    std::string ret;
    do
    {
        if (!ret.empty())
        {
            ret.insert(0, "/");
        }
        ret.insert(0, currentNode->name ? reinterpret_cast<const char*>(currentNode->name) : "#nullname");
        currentNode = currentNode->parent;
    }
    while(currentNode);

    return ret;
}

std::string XmlReader::getCurrentDebugContext() const
{
    xmlNodePtr currentNode = xmlTextReaderCurrentNode(xmlReader_);
    if (!currentNode)
    {
        return " no context";
    }

    const std::string ret = std::string(" line: ") + boost::lexical_cast<std::string>(currentNode->line) + " path:" + getXmlNodePath(currentNode);

    return ret;
}

const char *XmlReader::getName() const
{
    const char *ret = reinterpret_cast<const char*>(xmlTextReaderConstName(xmlReader_));
    if (!ret)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderConstName returned 0"));
    }
//    ISAAC_THREAD_CERR << "Name: " << ret << std::endl;
    return ret;
}

bool XmlReader::assertName(const char *name, const bool noThrow/* = false*/) const
{
    const char *currentName = getName();
    if (strcmp(name, currentName))
    {
        if (noThrow)
        {
            return false;
        }

        BOOST_THROW_EXCEPTION(xml::XmlReaderException(std::string("Unexpected element ") + currentName +
            " while looking for element " + name + getCurrentDebugContext()));
    }
    return true;
}

bool XmlReader::checkName(const char *name) const
{
    const char *currentName = getName();
    if (strcmp(name, currentName))
    {
        return false;
    }
    return true;
}

xmlReaderTypes XmlReader::getNodeType() const
{
    int nodeType = xmlTextReaderNodeType(xmlReader_);
    if (-1 == nodeType)
    {
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderNodeType failed"));
    }
    return static_cast<xmlReaderTypes>(nodeType);
}

bool XmlReader::read()
{
    int res = xmlTextReaderRead(xmlReader_);
    if (-1 == res)
    {
        good_ = false;
        BOOST_THROW_EXCEPTION(XmlReaderException("xmlTextReaderRead failed"));
    }
    ISAAC_ASSERT_MSG(0 == res || 1 == res, "Invalid result from xmlTextReaderRead : " << res);
    return 1 == res;
}


} // namespace xml
} // namespace isaac
