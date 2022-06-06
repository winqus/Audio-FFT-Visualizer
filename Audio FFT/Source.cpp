#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

/* Audio input and drawing */
#include <SFML/Graphics.hpp>
#include <SFML/System/Vector2.hpp>
#include <SFML/Window.hpp>
#include <SFML/Audio.hpp>

//#include <fftw3.h> // For FFT comparison and future use

using std::cout;
using std::endl;
using std::string;


template <typename T>
T clamp(const T& value, const T& minimum, const T& maximum) {
	return std::max(minimum, std::min(value, maximum));
}

class Slider : public sf::RectangleShape {
private:
	float m_statingValue;
	float m_value;
	float m_minValue;
	float m_maxValue;
	float m_step;
	bool m_isFocused;

	sf::Font m_font;
	sf::Text m_title, m_val;
	sf::RectangleShape m_path;
	sf::RectangleShape m_slider;

public:
	sf::Color sliderColor;
	bool showDecimal;

	Slider(sf::Vector2f size, float val, float minVal, float maxVal, float step, std::string title = "") : RectangleShape(size) {
		m_value = val;
		m_statingValue = m_value;
		m_minValue = std::min<float>(minVal, maxVal);
		m_maxValue = std::max<float>(minVal, maxVal);
		if (m_minValue == m_maxValue)
			m_minValue == 0.;
		m_step = step;
		m_isFocused = false;
		showDecimal = false;

		if (!m_font.loadFromFile("Lato-Black.ttf"))
			std::cout << "Couldn't load font for slider\n";


		this->setFillColor(sf::Color(31, 31, 31));
		this->setOutlineThickness(1.);
		this->setOutlineColor(sf::Color::White);

		setTitle(title);

		m_val.setCharacterSize(15);
		m_val.setFont(m_font);
		m_val.setFillColor(sf::Color::White);

		sliderColor = sf::Color(122, 125, 131);
	}
	void setTitle(std::string title) {
		m_title.setString(title);
		m_title.setCharacterSize(12);
		m_title.setFont(m_font);
		m_title.setFillColor(sf::Color::White);
		m_title.setOrigin(
			m_title.getGlobalBounds().width / 2.,
			m_title.getGlobalBounds().height
		);
		//m_title.setScale(0.5f, 0.5f);
	}
	void setValText() {
		m_val.setString(showDecimal ? std::to_string(m_value) : std::to_string((int)m_value));
		m_val.setOrigin(m_val.getGlobalBounds().width / 2., 0.f);
	}
	void SetValue(float value) {
		m_value = clamp<float>(value, m_minValue, m_maxValue);
	}
	void ChangeBySteps(int steps) {
		SetValue(m_value + m_step * steps);
	}
	void Update() {
		m_path = sf::RectangleShape(sf::Vector2f(this->getLocalBounds().width * 0.7f, this->getLocalBounds().height * 0.8f));
		m_path.setFillColor(sliderColor);
		m_path.setOrigin(m_path.getLocalBounds().width / 2., m_path.getLocalBounds().height / 2.);
		m_path.setPosition(
			this->getGlobalBounds().left + this->getLocalBounds().width / 2.,
			this->getGlobalBounds().top + this->getLocalBounds().height / 2.
		);
		m_title.setPosition(
			this->getGlobalBounds().left + this->getGlobalBounds().width / 2,
			this->getGlobalBounds().top - 8.f
		);
		setValText();
		m_val.setPosition(
			this->getGlobalBounds().left + this->getGlobalBounds().width / 2,
			this->getGlobalBounds().top + this->getGlobalBounds().height + 0.f
		);
		m_slider = sf::RectangleShape(m_path);
		m_slider.setFillColor(sf::Color(71, 201, 176));
		m_slider.setSize(sf::Vector2f(
			m_path.getSize().x,
			m_path.getSize().y * (m_value / (m_maxValue - m_minValue))
		));
		
	}
	bool UpdateFocus(sf::Vector2f mouseClickPos) {
		auto clickArea = sf::FloatRect(mouseClickPos.x, mouseClickPos.y, 0.1f, 0.1f);
		if (this->getGlobalBounds().intersects(clickArea)) {
			this->setOutlineColor(sf::Color::Red); // red
			if (m_isFocused && m_path.getGlobalBounds().intersects(clickArea)) {
				SetValue( 
					m_maxValue * abs(clickArea.top - m_path.getGlobalBounds().top) / m_path.getGlobalBounds().height
				);
			}
			m_isFocused = true;
		}
		else {
			m_isFocused = false;
			if(m_value == m_statingValue)
				this->setOutlineColor(sf::Color::White);
			else
				this->setOutlineColor(sf::Color(255, 233, 154)); // yellow
		}
		return m_isFocused;
	}
	bool CheckFocus() const {
		return m_isFocused;
	}
	float getValue() const {
		return m_value;
	}
	float getStartingValue() const {
		return m_statingValue;
	}
	void resetValue() {
		m_value = m_statingValue;
	}
	void draw(sf::RenderWindow& window) {
		Update();
		window.draw(*this);
		window.draw(m_path);
		window.draw(m_slider);
		window.draw(m_title);
		window.draw(m_val);
	}
};

class BarChart : public sf::RectangleShape {
private:
	std::vector<sf::RectangleShape> majorTickMarks, minorTickMarks;
	int majorTickSpacing, minorTickSpacing;
	sf::RectangleShape majorTickMarkTemplate, minorTickMarkTemplate;
	double m_scaleValue;
	sf::Font m_font;
	sf::Text m_title;
	sf::Text m_majorLabel;
public:
	std::vector<sf::RectangleShape> rects;
	sf::Vector2f barSizeMultiplierForOrigin;
	sf::Vector2f spacingVer, spacingHor;
	sf::Vector2f axisPlacementOffsetBySize_sizeMultiplier;
	sf::Color barColor;
	sf::Color borderColor;
	sf::Color tickColor;
	bool enableMajorTickMarks, enableMinorTickMarks;
	bool drawBorder;

	bool doAutoScale;
	bool doClipValues;

	bool useMajorTickLabels;
	bool majorTickLabelsDecimal;
	float majorTickLabelMultiplier;


	BarChart(float X, float Y, size_t dataSize, std::string title = "test title") {
		BarChart(sf::Vector2f(X, Y), dataSize, title);
	}
	BarChart(sf::Vector2f size, size_t dataSize, std::string title="test title") : RectangleShape(size) {
		setDefaults(dataSize);
		setTitle(title);
	}

	void setDefaults(size_t dataSize) {
		doClipValues = false;
		doAutoScale = false;
		drawBorder = true;
		useMajorTickLabels = true;
		majorTickLabelsDecimal = false;
		majorTickLabelMultiplier = 1.0;

		m_scaleValue = 1.0;

		if (dataSize < 10)
			dataSize = 10;
		barColor = sf::Color::White;
		borderColor = sf::Color(218U, 59U, 112U, 100U);
		tickColor = sf::Color::White;

		m_majorLabel = sf::Text("0", m_font, 8);
		m_majorLabel.setFillColor(sf::Color::White);
		majorTickSpacing = 32;
		minorTickSpacing = 8;


		this->setOrigin(this->getGlobalBounds().width / 2, 0.0f);
		this->setOutlineColor(borderColor);
		this->setOutlineThickness(3.0f);
		this->setFillColor(sf::Color(255U, 0U, 0U, 0U));

		if (!m_font.loadFromFile("Lato-Black.ttf"))
			std::cout << "Couldn't load font for bar chart\n";

		float placementWidth = this->getGlobalBounds().width - spacingHor.x - spacingHor.y;
		auto rectTemplate = sf::RectangleShape(sf::Vector2f(placementWidth / dataSize, 30.0f));
		rectTemplate.setOrigin(0.0f, rectTemplate.getGlobalBounds().height);

		spacingVer = sf::Vector2f(0.0f, 0.0f);
		spacingHor = spacingVer;

		axisPlacementOffsetBySize_sizeMultiplier = sf::Vector2f(0.0f, 0.0f);


		prepare(dataSize, rectTemplate, sf::Vector2f(0.0f, 1.0f));
		

	}

	void prepare(size_t dataSize, sf::RectangleShape rectTemplate, sf::Vector2f barSizeMultiplierForOrigin) {
		this->barSizeMultiplierForOrigin = barSizeMultiplierForOrigin;
		
		rects.resize(dataSize);
		for (size_t i = 0; i < rects.size(); i++) {
			auto& rect = rects[i];
			rect = rectTemplate;
		}

		setMajorTickMarks(true, majorTickSpacing, sf::RectangleShape(sf::Vector2f(1.0f, 16.0f)));
		setMinorTickMarks(true, minorTickSpacing, sf::RectangleShape(sf::Vector2f(1.0f, 8.0f)));
	}
	void plotData(const double data[], const int n) {
		int N = rects.size();
		if (N == 0 || n <= 0)
			return;
		float placementWidth = this->getGlobalBounds().width - spacingHor.x - spacingHor.y;
		float placementHeight = this->getSize().y - spacingVer.x - spacingVer.y;
		float minBarWidth = placementWidth / N;
		float barWidth = rects[0].getGlobalBounds().width >= minBarWidth ? rects[0].getGlobalBounds().width : minBarWidth;
		float verticalMult = placementHeight / this->getGlobalBounds().height;
		float horizontalMult = placementWidth / this->getGlobalBounds().width;

		int mult = 1;
		if (n > N)
			mult = n / N;

		double maxBarValue = 1.0;
		if (doAutoScale)
			for (size_t i = 0; i < N; i++)
				maxBarValue = data[i] > maxBarValue ? data[i] : maxBarValue;
			

		for (size_t i = 0; i < N; i++) {
			auto& bar = rects[i];
			float barValue = data[i*mult];
			if (i < n) {
				if (doAutoScale) 
					barValue = (barValue / maxBarValue);
				bar.setSize(sf::Vector2f(barWidth * horizontalMult, abs(barValue) * placementHeight * (axisPlacementOffsetBySize_sizeMultiplier.y) * verticalMult));
				
			}
			else bar.setSize(sf::Vector2f(barWidth, 0.0f));
			
			
			bar.setOrigin(
				bar.getSize().x * barSizeMultiplierForOrigin.x,
				barValue >= 0.f ? bar.getSize().y * barSizeMultiplierForOrigin.y : 0.0f
			);
			bar.setPosition(
				(this->getGlobalBounds().left + spacingHor.x + this->getSize().x * axisPlacementOffsetBySize_sizeMultiplier.x * horizontalMult) + barWidth * i * (n < N ? (float)N/n : 1.f),
				(this->getGlobalBounds().top + this->getOutlineThickness() + spacingVer.x + this->getSize().y * axisPlacementOffsetBySize_sizeMultiplier.y * verticalMult)
				//(this->getGlobalBounds().top + spacingVer.xk + this->getSize().y * axisPlacementOffsetBySize_sizeMultiplier.y)
			);

			if (enableMinorTickMarks) {
				minorTickMarks[i].setPosition(
					rects[i].getGlobalBounds().left + rects[i].getGlobalBounds().width / 2.0f,
					this->getGlobalBounds().top + this->getGlobalBounds().height
				);
			}
			if (enableMajorTickMarks) {
				majorTickMarks[i].setPosition(
					rects[i].getGlobalBounds().left + rects[i].getLocalBounds().width / 2.0f,
					this->getGlobalBounds().top + this->getLocalBounds().height
				);
			}
		}

		if (enableMajorTickMarks)
			majorTickMarks[majorTickMarks.size() - 1].setPosition(
				2.f * majorTickMarks[majorTickMarks.size()-2].getPosition().x - majorTickMarks[majorTickMarks.size()-3].getPosition().x,
				this->getGlobalBounds().top + this->getGlobalBounds().height
			);
		
		m_title.setPosition(
			this->getPosition()
		);
	}
	void setMajorTickMarks(bool enabled) {
		enableMajorTickMarks = enabled;
	}
	void setMajorTickMarks(bool enabled, int unit) {
		setMajorTickMarks(enabled, unit, majorTickMarkTemplate);
	}
	void setMajorTickMarks(bool enabled, int unit, sf::RectangleShape rectTemplate) {
		enableMajorTickMarks = enabled;
		majorTickSpacing = unit;
		majorTickMarkTemplate = rectTemplate;
		majorTickMarkTemplate.setOrigin(majorTickMarkTemplate.getSize().x / 2.0f, majorTickMarkTemplate.getSize().y * 0.0f);


		majorTickMarks.resize(rects.size()+1);
		for (int i = 0; i < majorTickMarks.size(); i++) {
			auto& tickMark = majorTickMarks[i];
			tickMark = majorTickMarkTemplate;
		}
	}
	void setMinorTickMarks(bool enabled) {
		enableMinorTickMarks = enabled;
	}
	void setMinorTickMarks(bool enabled, int unit) {
		setMinorTickMarks(enabled, unit, minorTickMarkTemplate);
	}
	void setMinorTickMarks(bool enabled, int unit, sf::RectangleShape rectTemplate) {
		enableMinorTickMarks = enabled;
		minorTickSpacing = unit;
		minorTickMarkTemplate = rectTemplate;
		minorTickMarkTemplate.setOrigin(minorTickMarkTemplate.getSize().x / 2.0f, minorTickMarkTemplate.getSize().y * 0.0f);

		minorTickMarks.resize(rects.size());
		for (size_t i = 0; i < minorTickMarks.size(); i++) {
			auto& tickMark = minorTickMarks[i];
			tickMark = minorTickMarkTemplate;
		}
		
	}
	void setTitle(std::string title) {
		m_title.setString(title);
		m_title.setCharacterSize(24);
		m_title.setFont(m_font);
		m_title.setFillColor(sf::Color::White);
		m_title.setOrigin(
			m_title.getGlobalBounds().width / 2.,
			m_title.getGlobalBounds().height + 15.f
		);
		//m_title.setScale(0.5f, 0.5f);
	}
	void setScaleValue(double value = 1.) {
		m_scaleValue = value;
	}
	void drawChart(sf::RenderWindow& window) {
		for (auto& rect : rects) {
			rect.setFillColor(barColor);
			window.draw(rect);
		}
			
		if (drawBorder) {
			this->setOutlineColor(borderColor);
			window.draw(*this);
		}


		if (enableMinorTickMarks) {
			for (size_t i = 0; i < minorTickMarks.size(); i += minorTickSpacing) {
				minorTickMarks[i].setFillColor(tickColor);
				window.draw(minorTickMarks[i]);
			}
		}
		if (enableMajorTickMarks) {
			
			for (size_t i = 0; i < majorTickMarks.size(); i += majorTickSpacing) {
				majorTickMarks[i].setFillColor(tickColor);
				window.draw(majorTickMarks[i]);

				if (!useMajorTickLabels)
					continue;
				m_majorLabel.setString( majorTickLabelsDecimal ? std::to_string(i * majorTickLabelMultiplier) : std::to_string(int(i*majorTickLabelMultiplier)));
				m_majorLabel.setOrigin(m_majorLabel.getGlobalBounds().width / 2., 0.f);
				m_majorLabel.setPosition(
					majorTickMarks[i].getPosition().x,
					majorTickMarks[i].getGlobalBounds().top + majorTickMarks[i].getGlobalBounds().height
				);
				window.draw(m_majorLabel);
			}
		}

		window.draw(m_title);

	}
};



sf::Mutex mutex;
class MyRecorder : public sf::SoundRecorder {

private:
	double* m_samples;
	size_t m_sampleSize;
	double m_minAmplitude;
	double m_maxAmplitude;
	double m_scaleValue;
	bool m_doAutoScale;

	virtual bool onStart()
	{
		// initialize whatever has to be done before the capture starts
		cout << ">>> Recording started with device: " << this->getDevice() << "\n";
		cout << "Proccessing interval is set to: " << (sf::Int32)std::round((m_sampleSize / (double)this->getSampleRate()) * 1000.)+1 << "ms\n";
		this->setProcessingInterval(sf::milliseconds(
			(sf::Int32)std::round( (m_sampleSize / (double)this->getSampleRate()) * 1000.) + 1
		));
		
		// return true to start the capture, or false to cancel it
		return true;
	}

	virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
	{
		sf::Lock lock(mutex);

		double average = 0.;
		double maxVal = 0.;
		for (size_t i = 0; i < m_sampleSize; i++) {
			if (i < sampleCount)
				m_samples[i] = (double)samples[i];
			else m_samples[i] = 0.;

			if (abs(m_samples[i]) < m_minAmplitude)
				m_samples[i] = 0.;
			if (abs(m_samples[i]) > m_maxAmplitude)
				m_samples[i] = clamp<double>(m_samples[i], -m_maxAmplitude, m_maxAmplitude);

			if (m_doAutoScale && abs(m_samples[i]) > maxVal)
				maxVal = abs(m_samples[i]);

			average += m_samples[i];
		}
		average /= std::min<double>(sampleCount, m_sampleSize);

		for (size_t i = 0; i < m_sampleSize; i++) {
			m_samples[i] -= average;
			if (m_doAutoScale)
				m_samples[i] /= maxVal;
			else
				m_samples[i] /= m_scaleValue;
			
		}

		// return true to continue the capture, or false to stop it
		return true;
	}

	virtual void onStop()
	{
		cout << ">>> Recording stopped <<<\n";
	}
public:
	MyRecorder(size_t bufferSize, bool autoScale = true) : SoundRecorder() {
		m_sampleSize = bufferSize;
		m_samples = new double[bufferSize];
		
		m_minAmplitude = 0.0;
		m_maxAmplitude = 1000000.0;
		m_scaleValue = 1. / 10000.;
		m_doAutoScale = autoScale;
	}
	~MyRecorder() {
		this->stop();
		delete[] m_samples;
		m_samples = nullptr;
	}
	size_t getDataSize() const {
		return m_sampleSize;
	}
	const double* getSamples() const {
		return m_samples;
	}
	void setAmplitudeRange(double minValue = 0.0, double maxValue = 1000000.0) {
		m_minAmplitude = minValue;
		m_maxAmplitude = maxValue;
	}
	void setScaleValue(double value=1.) {
		m_scaleValue = value;
	}
	void setAutoScale(bool state) {
		m_doAutoScale = state;
	}
};



void DFT(const int sampleCount, const int sampleRate, const double* samples, double* output) {

	/*** X1 x1 example ***/
		//std::complex<double> z1 = 1i;
		//cout << 0.707 * std::exp(-1i * 2. * M_PI * (1. * 1 / 8) ) << endl;
	/**********************/

	using namespace std::complex_literals;

	std::complex<double>* X = new std::complex<double>[sampleRate / 2 + 1]; \

	for (int k = 0; k < sampleRate / 2; k++) {
		//cout << k << "k\n";
		std::complex<double> Xk = 0.;
		for (int n = 0; n < sampleCount; n++) {
			//Xk += samples[n] * (std::exp(-1i * 2. * M_PI * ((double)k * n / sampleRate)));
			Xk += samples[n] * (std::exp(-1i * 2. * M_PI * ((double)k * n / sampleCount)));
		}
		//X[k] = ((Xk / (double)sampleCount));
		X[k] = Xk;
		// 
		//X[k] = sqrt( pow(std::real(X[k]), 2) + pow(std::imag(X[k]), 2));
		output[k] = std::abs(X[k]);
	}

	delete[] X;
}


unsigned ReverseBitOrder(unsigned int value, int bitCount) {
	unsigned newVal = 0;
	for (; bitCount > 0; bitCount--, value >>= 1) {
		newVal <<= 1;
		newVal |= (value & 0x1);
	}
	return newVal;
}

void BitReverseOrderData(std::complex<double>* input, std::complex<double>* output, int log2_size) {
	for (size_t i = 0; i < (1<<log2_size); i++)
		output[ReverseBitOrder(i, log2_size)] = input[i];
}

void FFT(std::complex<double>* samples, std::complex<double>* output, int log2_size) {
	using namespace std;

	int size = 1 << log2_size;
	complex<double>* X = new std::complex<double>[size];

	BitReverseOrderData(samples, X, log2_size);

	for (size_t s = 1; s <= log2_size; s++) {
		int N = 1 << s;

		complex<double> W = 1., Wm = exp(-1i * 2. * (M_PI / N));
		
		for (size_t i = 0; i < (N>>1); i++) {
			for (size_t k = i; k < size; k += N) {
				complex<double> wX_kN2 = W * X[k + (N>>1)];
				complex<double> X_k = X[k];

				X[k] = X_k + wX_kN2;
				X[k + (N>>1)] = X_k - wX_kN2;
			}
			W = W * Wm;
		}
	}

	for (size_t i = 0; i < size; i++) {
		output[i] = X[i];
		//output[i] = X[i] * (2. / size);
	}

	delete[] X;
}

void Hanning(std::complex<double>* data, int n) {
	for (int i = 0; i < n; i++) {
		data[i] *= 0.5 * (1 - cos(2 * M_PI * i / (n - 1)));
	}
}

int main() {
/* CONSTANTS */
	const int bufferSize = 1024<<0;
	const float chartSpacing = 70.0f; 
	const float sliderLeftOffset = 100.0f;
	const float sliderSpacing = 150.0f;
	const double maxDesiredFrequency = 5000;

/* Check if recording is available */
	if (!sf::SoundBufferRecorder::isAvailable())
	{
		cout << "audio capture is not available on this system" << endl;
		return 1;
	}
	
/*** SETUP RECORDING DEVICE ***/
	MyRecorder audioRecorder(bufferSize, true);

/* Choose input source */
	auto devices = audioRecorder.getAvailableDevices();
	cout << "Choose recording device:\n";
	for (size_t i = 1; i <= devices.size(); i++)
		cout << "\t" << i << ") " << devices[i - 1] << "\n";
	cout << "Default device (input): " << audioRecorder.getDefaultDevice() << "\n";
	int choice;
	cout << ">>> Enter choice (number): "; std::cin >> choice;
	choice = clamp<int>(--choice, 0, devices.size()-1);
	audioRecorder.setDevice(devices[choice]);
/* Begin recording */
	audioRecorder.start();
	double sampleRate = audioRecorder.getSampleRate();
	cout << "Sample rate: " << sampleRate << "\n";
	cout << "Buffer size: " << bufferSize << "\n";
	double freqRes = (double)audioRecorder.getSampleRate() / bufferSize;
	cout << "Frequency resolution: " << freqRes << "\n";

	

/*** RENDERING SETUP ***/
	sf::RenderWindow window(sf::VideoMode(1000, 1400), "Audio Visualizer");
	window.setFramerateLimit(30);

/* Chart for audio amplitude values */
	BarChart ampChart(sf::Vector2f(window.getSize().x - 2 * 20.0f, 300.0f), bufferSize, "Real Time Audio Samples");
	ampChart.axisPlacementOffsetBySize_sizeMultiplier = sf::Vector2f(0.0f, 0.5f);
	ampChart.setPosition(window.getSize().x / 2.f, 50.0f);
	ampChart.spacingVer = sf::Vector2f(10.0f, 10.0f);

/* Chart for frequency bin values */
	BarChart freqBinChart(sf::Vector2f(window.getSize().x - 2 * 20.0f, 300.0f), bufferSize, "Audio FFT Results");
	freqBinChart.barColor = sf::Color(230U, 145U, 56U);
	freqBinChart.axisPlacementOffsetBySize_sizeMultiplier = sf::Vector2f(0.0f, 1.0f);
	freqBinChart.setPosition(
		ampChart.getPosition().x, 
		ampChart.getGlobalBounds().top + ampChart.getGlobalBounds().height + chartSpacing
	);
	freqBinChart.spacingVer = sf::Vector2f(10.0f, 1.0f);
	freqBinChart.majorTickLabelMultiplier = freqRes;

/* Chart for frequency bin values */
	BarChart freqBinRangedChart(
		sf::Vector2f(window.getSize().x - 2 * 20.0f, 300.0f), int(maxDesiredFrequency / freqRes), "Frequencies 0-"+std::to_string(int(maxDesiredFrequency))+" Hz"
	);
	freqBinRangedChart.barColor = sf::Color(109, 173, 240);
	freqBinRangedChart.axisPlacementOffsetBySize_sizeMultiplier = sf::Vector2f(0.0f, 1.0f);
	freqBinRangedChart.setPosition(
		freqBinChart.getPosition().x,
		freqBinChart.getGlobalBounds().top + freqBinChart.getGlobalBounds().height + chartSpacing
	);
	freqBinRangedChart.setMajorTickMarks(true, 6);
	freqBinRangedChart.setMinorTickMarks(true, 3);
	freqBinRangedChart.majorTickLabelMultiplier = freqRes;

/* Sliders */
	Slider slidMinAmp(sf::Vector2f(20., 150.), 0., 0., 32000., 5., "MinAmp");
		slidMinAmp.setPosition(freqBinRangedChart.getGlobalBounds().left + sliderLeftOffset, freqBinRangedChart.getGlobalBounds().top + freqBinRangedChart.getLocalBounds().height + chartSpacing);
	Slider slidMaxAmp(sf::Vector2f(20., 150.), 50000., 1., 50000., 100., "MaxAmp");
	Slider slidScaleAmp(sf::Vector2f(20., 150.), 5000., 1., 50000., 100., "ScaleAmp (auto)");
	Slider slidScaleFreq(sf::Vector2f(20., 150.), 2., 0.5, 100., 0.5, "ScaleFreq");
		slidScaleFreq.showDecimal = true;
	Slider slidFreqSubVal(sf::Vector2f(20., 150.), 0., 0.0, 1., 0.01, "FreqSubVal");
		slidFreqSubVal.showDecimal = true;
	Slider slidFreqAutoScale(sf::Vector2f(20., 150.), 1., 0., 1., 1., "AutoScaleFreq");

	std::vector<Slider*> sliders = { 
		&slidMinAmp, &slidMaxAmp, &slidScaleAmp, &slidScaleFreq, &slidFreqSubVal, &slidFreqAutoScale
	};
	
	for (size_t i = 1; i < sliders.size(); i++)
		sliders[i]->setPosition(sliders[i - 1]->getGlobalBounds().left + sliderSpacing, freqBinRangedChart.getGlobalBounds().top + freqBinRangedChart.getLocalBounds().height + chartSpacing);

	std::complex<double> complexInput[bufferSize] = {};
	std::complex<double> complexOutput[bufferSize] = {};
	
	double freqBins[bufferSize] = {};
	double binsForRangedPlot[bufferSize/2] = {};
	
	//// For future use
	//fftw_plan plan = fftw_plan_dft_1d(bufferSize, reinterpret_cast<fftw_complex*>(&complexInput[0]), reinterpret_cast<fftw_complex*>(&complexOutput[0]), FFTW_FORWARD, FFTW_ESTIMATE);

/*** Application ***/
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::Resized)
			{
				// update the view to the new size of the window
				sf::FloatRect visibleArea(0.f, 0.f, event.size.width, event.size.height);
				window.setView(sf::View(visibleArea));

				ampChart.setPosition(window.getSize().x / 2.f, 50.0f);
				freqBinChart.setPosition(ampChart.getPosition().x, ampChart.getGlobalBounds().top + ampChart.getGlobalBounds().height + chartSpacing);
				freqBinRangedChart.setPosition(freqBinChart.getPosition().x, freqBinChart.getGlobalBounds().top + freqBinChart.getGlobalBounds().height + chartSpacing);
				
				slidMinAmp.setPosition(freqBinRangedChart.getGlobalBounds().left + sliderLeftOffset, freqBinRangedChart.getGlobalBounds().top + freqBinRangedChart.getLocalBounds().height + chartSpacing);
				for(size_t i = 1; i < sliders.size(); i++)
					sliders[i]->setPosition(sliders[i-1]->getGlobalBounds().left + sliderSpacing, freqBinRangedChart.getGlobalBounds().top + freqBinRangedChart.getLocalBounds().height + chartSpacing);
			}
			if (event.type == sf::Event::MouseButtonReleased) {
				if (event.mouseButton.button == sf::Mouse::Left) {
					for(auto slider : sliders)
						slider->UpdateFocus(window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y)));
					if (slidScaleAmp.CheckFocus()) {
						auto val = slidScaleAmp.getValue();
						if (val != slidScaleAmp.getStartingValue()) {
							audioRecorder.setAutoScale(false);
							audioRecorder.setScaleValue(val);
							slidScaleAmp.setTitle("ScaleAmp");
						}
					}
				}
				else if (event.mouseButton.button == sf::Mouse::Right) {
					for (auto slider : sliders)
						if(slider->UpdateFocus(window.mapPixelToCoords(sf::Vector2i(event.mouseButton.x, event.mouseButton.y))))
							slider->resetValue();

					if (slidScaleAmp.CheckFocus()) {
						audioRecorder.setAutoScale(true);
						audioRecorder.setScaleValue(slidScaleAmp.getValue());
						slidScaleAmp.setTitle("ScaleAmp (auto)");
					}
				}
				audioRecorder.setAmplitudeRange(slidMinAmp.getValue(), slidMaxAmp.getValue());
			}
			if (event.type == sf::Event::MouseWheelScrolled) {
				if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel) {
					bool slidersFocused = false;

					for(auto slider : sliders)
						if (slider->CheckFocus()) {
							slidersFocused = true;
							slider->ChangeBySteps(event.mouseWheelScroll.delta);
						}
					
					if (slidScaleAmp.CheckFocus()) {
						audioRecorder.setAutoScale(false);
						audioRecorder.setScaleValue(slidScaleAmp.getValue());
						slidScaleAmp.setTitle("ScaleAmp");
					}
					audioRecorder.setAmplitudeRange(slidMinAmp.getValue(), slidMaxAmp.getValue());
				}
			}
		}

		window.clear(sf::Color::Black);

		mutex.lock();

		auto audioSamples = audioRecorder.getSamples();
		for (size_t i = 0; i < bufferSize; i++) {
			complexInput[i] = audioSamples[i];
		}
		
		Hanning(complexInput, bufferSize); //// Apply window function to samples

		//fftw_execute(plan); //// faster, for future use (+uncomment plan and cleanup)
		FFT(complexInput, complexOutput, log2(bufferSize)); //// slower, for learning purposes

		double maxBin = 1.;
		for (size_t i = 0; i < bufferSize; i++) {
			freqBins[i] = abs(complexOutput[i]) * 2. / bufferSize; // Convert complex to pure real doubles

			if (slidFreqSubVal.getValue() > 0.) // Clamping
				freqBins[i] = clamp<double>(freqBins[i] - slidFreqSubVal.getValue(), 0., freqBins[i]);

			freqBins[i] *= slidScaleFreq.getValue(); // Scaling

			if(i < bufferSize / 2) // Only half is required (symmetrical data)
				maxBin = freqBins[i] > maxBin ? freqBins[i] : maxBin;

			if (i < int(maxDesiredFrequency / freqRes)) {
				binsForRangedPlot[i] = freqBins[i];
			}
		}
		/* Auto scalling frequency bins */
		if (slidFreqAutoScale.getValue() == 1.) {
			for (size_t i = 0; i < bufferSize; i++) {
				freqBins[i] /= maxBin;
				if (i < int(maxDesiredFrequency / freqRes)) {
					binsForRangedPlot[i] /= maxBin;
				}
			}
		}

		ampChart.plotData(audioRecorder.getSamples(), bufferSize);
		freqBinChart.plotData(freqBins, bufferSize);
		freqBinRangedChart.plotData(binsForRangedPlot, int(maxDesiredFrequency / freqRes));

		ampChart.drawChart(window);
		freqBinChart.drawChart(window);
		freqBinRangedChart.drawChart(window);

		for (auto slider : sliders)
			slider->draw(window);
		mutex.unlock();


		window.display();
	}
	//// For future use
	//fftw_destroy_plan(plan);
	//fftw_cleanup();
	return 0;
}