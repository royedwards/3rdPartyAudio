#ifndef EmbededWidgets_hxx
#define EmbededWidgets_hxx
#include <CLAM/Processing.hxx>
#include <QtGui/QWidget>
#include <map>

namespace CLAM
{

class EmbededWidgetCreatorBase
{
protected:
	typedef std::string Key;
private:
	typedef std::map<Key, EmbededWidgetCreatorBase*> Creators;
	static Creators & creators();
protected:
	EmbededWidgetCreatorBase(const Key & processingTypeName)
	{
		creators().insert(std::make_pair(processingTypeName, this));
	}
	virtual ~EmbededWidgetCreatorBase() {} // TODO: to the cxx
	virtual QWidget * concreteCreate(CLAM::Processing * processing, QWidget * parent) = 0;
public:
	static QWidget * create(CLAM::Processing * processing, QWidget * parent)
	{
		Key type = processing->GetClassName();
		Creators::iterator it = creators().find(type);
		if (it==creators().end()) return 0;
		return it->second->concreteCreate(processing, parent);
	}
};

template <typename WidgetType>
class EmbededWidgetCreator : public EmbededWidgetCreatorBase
{
public:
	EmbededWidgetCreator(const Key & type)
		: EmbededWidgetCreatorBase(type)
	{}
	virtual QWidget * concreteCreate(CLAM::Processing * processing, QWidget * parent)
	{
		return new WidgetType(processing);
	}
};

// Temporary duplication of classes with EmbededMonitorCreator as there are two 
// creation interfaces, passing the monitor on construction or in the setDataSource

template <typename WidgetType, typename MonitorType>
class EmbededMonitorCreator : public EmbededWidgetCreatorBase
{
public:
	EmbededMonitorCreator(const Key & type)
		: EmbededWidgetCreatorBase(type)
	{}
	virtual QWidget * concreteCreate(CLAM::Processing * processing, QWidget * parent)
	{
		MonitorType * monitor = dynamic_cast<MonitorType*>(processing);
		if (not monitor) return 0;
		WidgetType * widget = new WidgetType(parent);
		widget->setDataSource(*monitor);
		return widget;
	}
};

template <typename WidgetType, typename MonitorType>
class EmbededMonitorCreator2 : public EmbededWidgetCreatorBase
{
public:
	EmbededMonitorCreator2(const Key & type)
		: EmbededWidgetCreatorBase(type)
	{}
	virtual QWidget * concreteCreate(CLAM::Processing * processing, QWidget * parent)
	{
		MonitorType * monitor = dynamic_cast<MonitorType*>(processing);
		if (not monitor) return 0;
		return new WidgetType(parent,monitor);
	}
};

}

#endif//EmbededWidgets_hxx

