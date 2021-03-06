#include "InControlQtBinder.hxx"
#include "QtSlot2Control.hxx"
#include <QAbstractSlider>

static CLAM::InControlQtBinder inControlBinder;


CLAM::InControlQtBinder::InControlQtBinder()
	: _prefix("InControl__")
	, _controlProperty("clamInControl")
{
}

bool CLAM::InControlQtBinder::handles(QObject * uiElement)
{
	if (not isSubClassOf<QAbstractSlider>(uiElement)) return false;
	if (uiElement->objectName().startsWith(_prefix)) return true;
	if (not hasProperty(uiElement, _controlProperty)) return false;
	// TODO: Check the network control type
	return true;
}

static std::string networkName(QObject * uiElement, const QString & prefix, const char * property, QStringList & errors)
{
	if (uiElement->objectName().startsWith(prefix))
		return CLAM::QtBinder::widget2NetworkName(prefix,uiElement->objectName());

	if (CLAM::QtBinder::hasProperty(uiElement, property))
		return uiElement->property(property).toString().toStdString();

	CLAM::QtBinder::error(errors,
		QString("InControlQtBinder: "
			"Widget name '%1' should start by with '%2', "
			"or have the '%3' custom property defined")
			.arg(uiElement->objectName())
			.arg(prefix)
			.arg(property)
			);
	return "";
}

bool CLAM::InControlQtBinder::bind(QObject * slider, CLAM::Network & network, QStringList & errors)
{
	std::string controlName = networkName(
		slider, _prefix, _controlProperty, errors);

	if (controlName.empty()) return false;

	std::cout << "* Mapped widget using control bounds (map: 200:1->bounds): " << controlName << std::endl;

	if (missingInControl(controlName, network, errors)) return false;

	CLAM::InControlBase & receiver = network.GetInControlByCompleteName(controlName);

	if (receiver.GetTypeId() != typeid(CLAM::TControlData)) error(errors,
		QString("InControlQtBinder: '%1' is not a Float control.")
			.arg(controlName.c_str())
			);

	float widgetMinimum = propertyDefault( slider, "minimum", 0);
	float widgetMaximum = propertyDefault( slider, "maximum", 100);

	float minimum = propertyDefault( slider, "clamMinimum",
		receiver.IsBounded()? receiver.LowerBound() : widgetMinimum);
	float maximum = propertyDefault( slider, "clamMaximum",
		receiver.IsBounded()? receiver.UpperBound() : widgetMaximum);

	std::cout
		<< "Mapping bounds " << widgetMinimum << "," << widgetMaximum
		<< " to " << minimum << "," << maximum << std::endl;

	QtSlot2Control * notifier = new QtSlot2Control(
		controlName.c_str(), 
		minimum, maximum,
		widgetMinimum, widgetMaximum);
	notifier->setParent(slider);
	notifier->linkControl(receiver);
	notifier->connect(
		slider,SIGNAL(valueChanged(int)),
		SLOT(sendMappedControl(int)));
	return true;
}

void CLAM::InControlQtBinder::setPropertyIfPresent(QObject * uiElement, const char * property, int value)
{
	if (hasProperty(uiElement, property))
		uiElement->setProperty(property, QVariant(value));
}


class BoolInControlQtBinder : public CLAM::QtBinder
{
	QString _prefix;
public:
	BoolInControlQtBinder()
		: _prefix("InControlBool__")
	{
	}
	bool handles(QObject * uiElement)
	{
		// TODO: check that it has a togled signal
		if (not isSubClassOf<QWidget>(uiElement)) return false;
		if (uiElement->objectName().startsWith(_prefix)) return true;
		return false;
	}
	bool bind(QObject * uiElement, CLAM::Network & network, QStringList & errors)
	{
		std::string controlName=widget2NetworkName(_prefix,uiElement->objectName());
		std::cout << "* Bool Control (to float map): " << controlName << std::endl;

		if (missingInControl(controlName, network, errors)) return false;
		CLAM::InControlBase & receiver = network.GetInControlByCompleteName(controlName);
		QtSlot2Control * notifier = new QtSlot2Control(controlName.c_str());
		notifier->setParent(uiElement);
		notifier->linkControl(receiver);
		notifier->connect(
				uiElement, SIGNAL(toggled(bool)),
				SLOT(sendBooleanControl(bool)));
		return true;
	}
};

static BoolInControlQtBinder boolBinder;

class IntegerInControlQtBinder : public CLAM::QtBinder
{
	QString _prefix;
public:
	IntegerInControlQtBinder()
		: _prefix("InControlInteger__")
	{
	}
	bool handles(QObject * uiElement)
	{
		// TODO: check that it has a togled signal
		if (not isSubClassOf<QWidget>(uiElement)) return false;
		if (uiElement->objectName().startsWith(_prefix)) return true;
		return false;
	}
	bool bind(QObject * uiElement, CLAM::Network & network, QStringList & errors)
	{
		std::string controlName=widget2NetworkName(_prefix,uiElement->objectName());
		std::cout << "* Integer Control (to float map): " << controlName << std::endl;

		if (missingInControl(controlName, network, errors)) return false;
		CLAM::InControlBase & receiver = network.GetInControlByCompleteName(controlName);
		QtSlot2Control * notifier = new QtSlot2Control(controlName.c_str());
		notifier->setParent(uiElement);
		notifier->linkControl(receiver);
		notifier->connect(
				uiElement, SIGNAL(valueChanged(int)),
				SLOT(sendControl(int)));
		return true;
	}
};

static IntegerInControlQtBinder integerBinder;

class FloatInControlQtBinder : public CLAM::QtBinder
{
	QString _prefix;
public:
	FloatInControlQtBinder()
		: _prefix("InControlFloatr__")
	{
	}
	bool handles(QObject * uiElement)
	{
		// TODO: check that it has a togled signal
		if (not isSubClassOf<QWidget>(uiElement)) return false;
		if (uiElement->objectName().startsWith(_prefix)) return true;
		return false;
	}
	bool bind(QObject * uiElement, CLAM::Network & network, QStringList & errors)
	{
		std::string controlName =QtBinder::widget2NetworkName("InControlFloat__",uiElement->objectName());
		std::cout << "* Mapped Control (100:1): " << controlName << std::endl;

		if (missingInControl(controlName, network, errors)) return false;
		CLAM::InControlBase & receiver = network.GetInControlByCompleteName(controlName);
		QtSlot2Control * notifier = new QtSlot2Control(controlName.c_str());
		notifier->setParent(uiElement);
		notifier->linkControl(receiver);
		notifier->connect(
				uiElement,SIGNAL(valueChanged(int)),
				SLOT(sendMappedControl(int)));
		return true;
	}
};

static FloatInControlQtBinder floatBinder;

#include <CLAM/ControlSource.hxx>
#include <QtGui/QAbstractSlider>

void ControlSourceSender::send(int value)
{
	_source->Do(value/100.);
}

namespace CLAM { class ControlSourceQtBinder; }

class CLAM::ControlSourceQtBinder : public CLAM::QtBinder
{
	std::list<ControlSourceSender *> _toDelete;
	QString _prefix;
public:
	ControlSourceQtBinder()
		: _prefix("ControlSource__")
	{
	}
	~ControlSourceQtBinder()
	{
		for (std::list<ControlSourceSender *>::iterator it=_toDelete.begin();
			it!=_toDelete.end(); it++)
			delete *it;
	}
	bool handles(QObject * uiElement)
	{
		// TODO: check that it has a togled signal
		if (not isSubClassOf<QAbstractSlider>(uiElement)) return false;
		if (uiElement->objectName().startsWith(_prefix)) return true;
		return false;
	}
	bool bind(QObject * uiElement, CLAM::Network & network, QStringList & errors)
	{
		std::string processingName=widget2NetworkName(_prefix,uiElement->objectName());
		std::cout << "* Slider connected to ControlSource port " << processingName << std::endl;
		Processing & processing = network.GetProcessing(processingName);
		ControlSource * source = dynamic_cast<ControlSource*> (&processing);
		if (not source) return error(errors,
			QString("ControlSourceQtBinder: Processing '%1' is of class '%2' "
				"And not 'ControlSource' as expected")
				.arg(processingName.c_str())
				.arg(processing.GetClassName())
				);
		ControlSourceSender * bind = new ControlSourceSender(source);
		bind->setParent(uiElement);
		_toDelete.push_back(bind);
		QObject::connect(
			uiElement, SIGNAL(valueChanged(int)),
			bind, SLOT(send(int)));
		return true;
	}
};

static CLAM::ControlSourceQtBinder controlSourceBinder;


