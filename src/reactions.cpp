#include "reactions.h"

#include "spark/collisions/reactions/he_reactions.h"
#include "rapidcsv.h"

namespace {
spark::collisions::CrossSection load_cross_section(const std::filesystem::path& path, double energy_threshold) {
    spark::collisions::CrossSection cs;
    rapidcsv::Document doc(path.string(), rapidcsv::LabelParams(-1, -1),
                           rapidcsv::SeparatorParams(';'));
    cs.energy = doc.GetColumn<double>(0);
    cs.cross_section = doc.GetColumn<double>(1);
    cs.threshold = energy_threshold;
    return cs;
}
}  // namespace

spark::collisions::Reactions<2, 3> spark::reactions::load_electron_reactions(
    const std::filesystem::path& dir,
    const Parameters& par,
    spark::particle::ChargedSpecies<2, 3>& ions) {
    spark::collisions::Reactions<2, 3> electron_reactions;
    electron_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeElectronElasticCollision<2, 3>>(
            spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Elastic_He.csv", 0.0)));

    electron_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeExcitationCollision<2, 3>>(
            spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Excitation1_He.csv", 19.82)));

    electron_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeExcitationCollision<2, 3>>(
            spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Excitation2_He.csv", 20.61)));

    electron_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeIonizationCollision<2, 3>>(
            ions, par.tg, spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Ionization_He.csv", 24.59)));

    return electron_reactions;
}

spark::collisions::Reactions<2, 3> spark::reactions::load_ion_reactions(const std::filesystem::path& dir,
                                                                   const Parameters& par) {
    spark::collisions::Reactions<2, 3> ion_reactions;
    ion_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeIonElasticCollision<2, 3>>(
            spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Isotropic_He.csv", 0.0)));

    ion_reactions.push_back(
        std::make_unique<spark::collisions::reactions::HeIonChargeExchangeCollision<2, 3>>(
            spark::collisions::reactions::HeCollisionConfig{par.m_he},
            load_cross_section(dir / "Backscattering_He.csv", 0.0)));

    return ion_reactions;
}