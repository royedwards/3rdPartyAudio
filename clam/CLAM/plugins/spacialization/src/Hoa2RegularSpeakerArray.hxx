#ifndef Hoa2RegularSpeakerArray_hxx
#define Hoa2RegularSpeakerArray_hxx
#include <CLAM/AudioInPort.hxx>
#include <CLAM/AudioOutPort.hxx>
#include <CLAM/Processing.hxx>
#include <CLAM/Audio.hxx>
#include <CLAM/InControl.hxx>
#include <CLAM/Filename.hxx>
#include <CLAM/Enum.hxx>
#include "Orientation.hxx"
#include "SpeakerLayout.hxx"
#include <cmath>

/**
 Decodes a Nth order signal for a regular array of speakers.
 Regular arrays can be decoded with a single constant for all the components 
 of the same degree, given an order of ambisonics.
 @param Order [Config] The order of ambisonics that it will be fed with.
 @param SpeakerLayout [Config] A file containing the target speaker layout.
 @param[out] XX [Port] Ambisonics component were XX is one of W,X,Y,Z,R,S,T,U,V,K,L,M,N,O,P,Q.
 @param[out] AAA [Ports] Audio signals to be emitted by the speaker AAA. \
   Here AAA is the label in the layout file or 01, 02... if IgnoreLabels is set (the default)
 @see SpeakerLayout for a description of the file format.
 @todo Make the order configurable
 @todo Make the decoding configurable
 @ingroup SpatialAudio
 @see AmbisonicsConventions
*/
class Hoa2RegularSpeakerArray : public CLAM::Processing
{
public:
	class DecodingCriteria : public CLAM::Enum
	{
	public:
		DecodingCriteria() : CLAM::Enum(ValueTable(), eInPhase2D) {}
		DecodingCriteria(tValue v) : CLAM::Enum(ValueTable(), v) {}
		DecodingCriteria(const std::string & s) : CLAM::Enum(ValueTable(), s) {}
		virtual CLAM::Component * Species() const {return new DecodingCriteria; }
		typedef enum {
			eInPhase2D=0,
			eInPhase3D=1,
		} tEnum;
		static tEnumValue * ValueTable()
		{
			static tEnumValue sValueTable[] =
			{
				{eInPhase2D,"In-phase 2D"},
				{eInPhase3D,"In-phase 3D"},
				{0,NULL}
			};
			return sValueTable;
		}
	};
	class Config : public CLAM::ProcessingConfig
	{
		DYNAMIC_TYPE_USING_INTERFACE( Config, 4, ProcessingConfig );
		DYN_ATTRIBUTE( 0, public, unsigned, Order);
		DYN_ATTRIBUTE( 1, public, CLAM::InFilename, SpeakerLayout);
		DYN_ATTRIBUTE( 2, public, bool, IgnoreLabels);
		DYN_ATTRIBUTE( 3, public, DecodingCriteria, DecodingCriteria);
	protected:
		void DefaultInit()
		{
			AddAll();
			UpdateData();
			SetOrder(1);
			SetIgnoreLabels(true);
		}
	};
private:
	SpeakerLayout _layout;
	typedef std::vector<CLAM::AudioOutPort*> OutPorts;
	typedef std::vector<CLAM::AudioInPort*> InPorts;
	OutPorts _outputs;
	InPorts _inputs;
	Config _config;
	enum {MaxSupportedOrder=3};
	double _decoding[MaxSupportedOrder+1];

public:
	const char* GetClassName() const { return "Hoa2RegularSpeakerArray"; }
	Hoa2RegularSpeakerArray(const Config& config = Config()) 
	{
		Configure( config );
	}
	const CLAM::ProcessingConfig & GetConfig() const
	{
		return _config;
	}

	bool ConcreteConfigure(const CLAM::ProcessingConfig& config)
	{
		CopyAsConcreteConfig(_config, config);
		unsigned buffersize = BackendBufferSize();
		unsigned order = _config.GetOrder();
		if (not _config.HasDecodingCriteria())
		{
			_config.AddDecodingCriteria();
			_config.UpdateData();
			_config.SetDecodingCriteria(DecodingCriteria::eInPhase2D);
		}
		bool ok = true;
		if (order>MaxSupportedOrder)
		{
			// If never configured just initialize a sane default to keep connections
			if (_inputs.size() == 0) ResizePortsToOrder(1, buffersize);
			ok = AddConfigErrorMessage("Ambisonics orders beyond 3rd are not supported");
			// Don't exit yet, we want to keep outports connections when loading from a network
		}
		else
		{
			ResizePortsToOrder(order, buffersize);
			ComputeDecoding(order);
		}

		std::string errorMessage;
		if (not _layout.load(_config.GetSpeakerLayout(), errorMessage))
			return AddConfigErrorMessage(errorMessage);
		else ResizePortsToLayout(buffersize);

		return ok;
	}
	bool Do()
	{
		// Ambisonics definition assures us one component for order 0 so we can do that:
		const unsigned nSamples = _inputs[0]->GetAudio().GetBuffer().Size();
		const unsigned nComponents = _inputs.size();
		const unsigned nSpeakers = _outputs.size();
		const CLAM::TData* components[nComponents];
		for (unsigned component=0; component<nComponents; component++)
			components[component] = &_inputs[component]->GetAudio().GetBuffer()[0];
		CLAM::TData* speakers[nSpeakers];
		for (unsigned speaker=0; speaker<nSpeakers; speaker++)
		{
			CLAM::Audio & audio = _outputs[speaker]->GetAudio();
			speakers[speaker] = &audio.GetBuffer()[0];
		}

		double componentWeight[nComponents];
		CLAM::SphericalHarmonicsDefinition * sh = CLAM::Orientation::sphericalHarmonics();
		for (unsigned speaker=0; speaker<nSpeakers; speaker++)
		{
			for (unsigned hoaComponent=0; hoaComponent<nComponents; hoaComponent++)
				componentWeight[hoaComponent] = _decoding[sh[hoaComponent].zProjection] * _layout.orientation(speaker).sphericalHarmonic(sh[hoaComponent]) / nSpeakers;
			CLAM::TData * speakerBuffer = speakers[speaker];
			for (unsigned sample=0; sample<nSamples; sample++)
			{
				double sampleValue = 0;
				for (unsigned component=0; component<nComponents; component++)
					sampleValue += components[component][sample] * componentWeight[component];
				speakerBuffer[sample] = 2 * sampleValue; //The 2 factor makes the volume comparable to VBAP
			}
		}
		for (unsigned speaker=0; speaker<nSpeakers; speaker++)
			_outputs[speaker]->Produce();
		for (unsigned component=0; component<nComponents; component++)
			_inputs[component]->Consume();
		return true;
	}
	~Hoa2RegularSpeakerArray()
	{
		for (unsigned speaker=0; speaker<_outputs.size(); speaker++)
			delete _outputs[speaker];
		for (unsigned component=0; component<_inputs.size(); component++)
			delete _inputs[component];
	}
private:
	void ResizePortsToOrder(unsigned order, unsigned buffersize)
	{
		// Set up the inputs according to ambisonics order
		CLAM::SphericalHarmonicsDefinition *sh = CLAM::Orientation::sphericalHarmonics();
		unsigned i=0;
		for (;sh[i].name; i++)
		{
			if (sh[i].order > order) break;
			if (i<_inputs.size())
			{
				_inputs[i]->SetSize( buffersize );
				_inputs[i]->SetHop( buffersize );
				continue;
			}
			CLAM::AudioInPort * port = new CLAM::AudioInPort( sh[i].name, this);
			port->SetSize( buffersize );
			port->SetHop( buffersize );
			_inputs.push_back( port );
		}
		unsigned actualSize=i;
		for (;i<_inputs.size(); i++)
			delete _inputs[i];
		_inputs.resize(actualSize);
	}
	void ResizePortsToLayout(unsigned buffersize)
	{
		// Set up the outputs according to the layout
		unsigned speakerToUpdate = firstDirtySpeaker();
		// delete existing speakers from the first one with different name
		for (unsigned speaker=0; speaker<speakerToUpdate; speaker++)
		{
			// Update the size and hop just in case
			_outputs[speaker]->SetSize( buffersize );
			_outputs[speaker]->SetHop( buffersize );
		}
		for (unsigned speaker=speakerToUpdate ; speaker<_outputs.size(); speaker++)
			delete _outputs[speaker];
		_outputs.resize(speakerToUpdate);
		// adding new speakers
		for (unsigned speaker=speakerToUpdate; speaker<_layout.size(); speaker++)
		{
			CLAM::AudioOutPort * port = new CLAM::AudioOutPort( portName(speaker), this);
			port->SetSize( buffersize );
			port->SetHop( buffersize );
			_outputs.push_back( port );
		}
	}
	std::string portName(unsigned speaker) const
	{
		if (_config.HasIgnoreLabels() and _config.GetIgnoreLabels())
		{
			std::ostringstream os;
			os << std::setw(2) << std::setfill('0') << (speaker+1);
			return os.str();
		}
		return _layout.name(speaker);
	}
	unsigned firstDirtySpeaker() const
	{
		for (unsigned speaker = 0; speaker<_layout.size(); speaker++)
		{
			if (speaker>=_outputs.size()) return speaker;
			if (_outputs[speaker]->GetName() != portName(speaker)) return speaker;
		}
		return _layout.size();
	}
	void ComputeDecoding(unsigned order)
	{
		DecodingCriteria criteria = _config.GetDecodingCriteria();
		double (*decodingFunction)(unsigned, unsigned) =
			(criteria==DecodingCriteria::eInPhase2D ? &inphaseDecoding2D :
			(criteria==DecodingCriteria::eInPhase3D ? &inphaseDecoding3D :
			/* default */ &inphaseDecoding2D ));

		for (unsigned i=0; i<=order; i++)
			_decoding[i] = decodingFunction(order, i);
		for (unsigned i=order+1; i<= MaxSupportedOrder; i++)
			_decoding[i] = 0;
	}
	static double inphaseDecoding2D(unsigned maxOrder, unsigned order)
	{
		/*
			1..n * 1..n / (1..n-l) / (1..n+l) = n-l+1..n / n+1..n+l
		*/
		double g = 1;
		for (unsigned i=maxOrder-order+1; i<=maxOrder; i++) g *= i;
		for (unsigned i=maxOrder+1; i<=maxOrder+order; i++) g /= i;
		if (order) g *= 2;
		return g;
	}
	static double inphaseDecoding3D(unsigned maxOrder, unsigned order)
	{
		double g=1;
		// TODO: Optimize to two productories as 2D decoding but be careful with spacial cases
		for (unsigned i=1; i<=maxOrder; i++) g*=i;
		for (unsigned i=1; i<=maxOrder+1; i++) g*=i;
		for (unsigned i=1; i<=maxOrder+order+1; i++) g/=i;
		for (unsigned i=1; i<=maxOrder-order; i++) g/=i;
		return g;
	}

};
#endif
