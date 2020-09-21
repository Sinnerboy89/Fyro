#include <JuceHeader.h>

int main(int argc, char* argv[])
{

    // parameters
    const float c = 343; // speed of sound in metres per second
    const float theta = juce::MathConstants<float>::pi / 4.0f; // azimuth of source
    const float gamma_le = -juce::MathConstants<float>::pi / 2.0f; // azimuth of left ear
    const float gamma_re = juce::MathConstants<float>::pi / 2.0f; // azimuth of right ear
    const float a_h = 0.0875f; // average radius of adult human head (according to B & D)
    const float e2e = 0.215f; // distance between ears, used for ITDs
    const float alpha_min = 0.1f; // values chosen by B& D as reasonable starting points before personalization
    const float theta_min = 150.0f * juce::MathConstants<float>::pi / 180.0f; // values chosen by B& D as reasonable starting points before personalization

    // load audio file (convert to mono if required)
    float srate = 0.0f;
    juce::AudioSampleBuffer fileBuffer;
    juce::AudioFormatManager formatManager;
    juce::FileChooser chooser("Select an audio file", {}, "*.wav");
    if (chooser.browseForFileToOpen()) {
        auto file = chooser.getResult();
        std::unique_ptr<juce::AudioFormatReader> reader(formatManager.createReaderFor(file));
        if (reader.get() != nullptr) {
            if (reader->numChannels == 2) {
                // TODO: mixdown stereo to mono
            }
            srate = (float)reader->sampleRate;
            fileBuffer.setSize(1, (int)reader->lengthInSamples);
            reader->read(&fileBuffer, 0, (int)reader->lengthInSamples, 0, true, true);
        }
    }

    // Calculate AOIs and ITDs
    const float w0 = c / a_h;
    float theta_l, theta_r, Td_l, Td_r;
    if (theta < 0) { // i.e. source is in front - left quarter space
        theta_l = -gamma_le + theta; // the left ear is the ipsilateral ear
        theta_r = gamma_re - theta; // the right ear is the contralateral ear
        Td_l = -(e2e / (2 * c)) * juce::dsp::FastMathApproximations::cos(theta_l) + (e2e / (2 * c)); // delay between centre and left ear
        Td_r = (e2e / (2 * c)) * (abs(theta_r) - (juce::MathConstants<float>::pi / 2)) + (e2e / (2 * c)); // delay between centre and right ear
    }
    else { // i.e.source is in front - right quarter space
        theta_l = -gamma_le + theta; // vice versa
        theta_r = gamma_re - theta; // vice versa
        Td_l = (e2e / (2 * c)) * (abs(theta_l) - (juce::MathConstants<float>::pi / 2)) + (e2e / (2 * c)); // time taken for signal to reach left ear
        Td_r = -(e2e / (2 * c)) * juce::dsp::FastMathApproximations::cos(theta_r) + (e2e / (2 * c)); // time taken for signal to reach right ear
    }

    // Calculate alphas for each ear
    const float alpha_l = (1 + alpha_min / 2) + (1 - alpha_min / 2) * juce::dsp::FastMathApproximations::cos(theta_l * juce::MathConstants<float>::pi / theta_min);
    const float alpha_r = (1 + alpha_min / 2) + (1 - alpha_min / 2) * juce::dsp::FastMathApproximations::cos(theta_r * juce::MathConstants<float>::pi / theta_min);

    // Set up ILD IIRs for each ear
    auto ild_coeffs_l = juce::dsp::IIR::Coefficients<float>::Coefficients(w0 + (alpha_l * srate), w0 - (alpha_l * srate), w0 + srate, w0 - srate);
    auto ild_coeffs_r = juce::dsp::IIR::Coefficients<float>::Coefficients(w0 + (alpha_r * srate), w0 - (alpha_r * srate), w0 + srate, w0 - srate);

    // Process audio

    // Save audio file


    return 0;
}
