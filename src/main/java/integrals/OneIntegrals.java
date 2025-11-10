package integrals;

import molecule.Atom;
import molecule.MoleculeIntegral;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import parameter.Parameter;

import java.util.List;

public class OneIntegrals {
    private static final Logger log = LoggerFactory.getLogger(OneIntegrals.class);

    public record AtomData(
            double orbC1, double orbC2,
            Atom.Coordinates coord1, Atom.Coordinates coord2, Atom.Coordinates coordNuc,
            Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2
    ) {}

    public record AtomCoordData(
            double orbC1, double orbC2,
            double coord1, double coord2, double nucCoord
    ) {}

    private final double[][] overlapInt;
    private final double[][] kineticInt;
    private final double[][] potentialInt;

    // Constructor — builds overlap integral matrix automatically
    public OneIntegrals(MoleculeIntegral molInt) {
        log.info("Overlap integrals computed");
        this.overlapInt = computeOverlap(molInt);

        log.info("Kinetic integrals computed");
        this.kineticInt = computeKinetic(molInt);

        log.info("Potential integrals computed");
        this.potentialInt = computePotential(molInt);
    }

    // Getter integrals matrix
    public double[][] getOverlapInt() {
        return overlapInt;
    }
    public double[][] getKineticInt() {
        return kineticInt;
    }
    public double[][] getPotentialInt() {
        return potentialInt;
    }

    // Compute overlap integrals matrix
    public static double[][] computeOverlap(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intOvelap = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intOvelap[i][j] = overlap(atomInt.get(i),atomInt.get(j));
                intOvelap[j][i] = intOvelap[i][j];
            }
        }
        return intOvelap;
    }

    private static double overlap(Atom.AtomIntegrals atomInt1,
                                  Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * ovInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    public static double ovInt(double orbC1, double orbC2,
                               Atom.Coordinates coord1, Atom.Coordinates coord2,
                               Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {

        double producto = 1.0;

        // Product over x, y, z components
        producto *= sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x());
        producto *= sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y());
        producto *= sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        // Compute distance squared between coord1 and coord2
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double eAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return eAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * producto;
    }

    private static double sRecurr(double alpha1, double beta1,
                                  double coord1, double coord2,
                                  int angCoord1coord, int angCoord2coord) {

        // Base cases
        if (angCoord1coord == 0 && angCoord2coord == 0) {
            return 1.0;
        } else if (angCoord1coord == 1 && angCoord2coord == 0) {
            return -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)));
        } else if (angCoord1coord > 1 && angCoord2coord == 0) {
            double term1 = -(coord1 - ((alpha1 * coord1 + beta1 * coord2) / (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 1, angCoord2coord);
            double term2 = ((angCoord1coord - 1) / (2.0 * (alpha1 + beta1)))
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord - 2, angCoord2coord);
            return term1 + term2;
        } else {
            return sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1)
                    + (coord1 - coord2)
                    * sRecurr(alpha1, beta1, coord1, coord2, angCoord1coord, angCoord2coord - 1);
        }
    }

    // Compute kinetic integrals matrix
    public static double[][] computeKinetic(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        int nOrb = atomInt.size();

        double[][] intKinetic = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intKinetic[i][j] = kinetic(atomInt.get(i),atomInt.get(j));
                intKinetic[j][i] = intKinetic[i][j];
            }
        }
        return intKinetic;
    }

    private static double kinetic(Atom.AtomIntegrals atomInt1,
                                  Atom.AtomIntegrals atomInt2) {

        double resultado = 0.0;

        for (int i = 0; i < atomInt1.exponent().size(); i++) {
            for (int j = 0; j < atomInt2.exponent().size(); j++) {
                resultado += normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                        * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                        * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                        * kinInt(atomInt1.exponent().get(i), atomInt2.exponent().get(j), atomInt1.coord(),
                        atomInt2.coord(), atomInt1.cartAngular(), atomInt2.cartAngular());
            }
        }

        return resultado;
    }

    private static double kinInt(double orbC1, double orbC2,
                                 Atom.Coordinates coord1, Atom.Coordinates coord2,
                                 Atom.AngCoordinates angCoord1, Atom.AngCoordinates angCoord2) {


        double producto1 = kiRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto2 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * kiRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * sRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double producto3 = sRecurr(orbC1, orbC2, coord1.x(), coord2.x(), angCoord1.x(), angCoord2.x())
                * sRecurr(orbC1, orbC2, coord1.y(), coord2.y(), angCoord1.y(), angCoord2.y())
                * kiRecurr(orbC1, orbC2, coord1.z(), coord2.z(), angCoord1.z(), angCoord2.z());

        double suma = producto1 + producto2 + producto3;

        // Distance squared between centers
        double dist2 = Math.pow(coord1.x() - coord2.x(), 2) +
                Math.pow(coord1.y() - coord2.y(), 2) +
                Math.pow(coord1.z() - coord2.z(), 2);

        double energyAB = Math.exp(-(orbC1 * orbC2 / (orbC1 + orbC2)) * dist2);

        return energyAB * Math.pow(Math.PI / (orbC1 + orbC2), 1.5) * suma;
    }

    private static double kiRecurr(double alpha, double beta,
                                   double coord1, double coord2,
                                   int angCoord1coord, int angCoord2coord) {

        double resultado;

        if (angCoord1coord == 0 && angCoord2coord == 0) {
            resultado = 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, 1);
        } else if (angCoord1coord > 0 && angCoord2coord == 0) {
            resultado = -angCoord1coord * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, 1);
        } else if (angCoord1coord == 0 && angCoord2coord > 0) {
            resultado = -angCoord2coord * alpha * sRecurr(alpha, beta, coord1, coord2, 1, angCoord2coord - 1)
                    + 2.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, 1, angCoord2coord + 1);
        } else {
            resultado = (
                    (angCoord1coord * angCoord2coord * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, angCoord2coord - 1))
                            - (2.0 * angCoord1coord * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord - 1, angCoord2coord + 1))
                            - (2.0 * angCoord2coord * alpha * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, angCoord2coord - 1))
                            + (4.0 * alpha * beta * sRecurr(alpha, beta, coord1, coord2, angCoord1coord + 1, angCoord2coord + 1))
            ) / 2.0;
        }

        return resultado;
    }

    // Compute potential integrals matrix
    public static double[][] computePotential(MoleculeIntegral molInt) {

        List<Atom.AtomIntegrals> atomInt = molInt.getAtoms();
        List<Integer> numAtomic = molInt.getAtNumbers();
        List<Atom.Coordinates> nucCoord = molInt.getNucCoord();
        int nOrb = atomInt.size();

        double[][] intPotential = new double[nOrb][nOrb];

        for (int i = 0; i < nOrb; i++) {
            for (int j = i; j < nOrb; j++) {
                intPotential[i][j] = potential(atomInt.get(i), atomInt.get(j), numAtomic, nucCoord);
                intPotential[j][i] = intPotential[i][j];
            }
        }
        return intPotential;
    }

    private static double potential(Atom.AtomIntegrals atomInt1,
                                  Atom.AtomIntegrals atomInt2, List<Integer> numAtomic, List<Atom.Coordinates> nucCoord) {

        double resultado = 0.0;
        for (int n = 0; n < nucCoord.size(); n++) {
            for (int i = 0; i < atomInt1.exponent().size(); i++) {
                for (int j = 0; j < atomInt2.exponent().size(); j++) {
                    AtomData atData = new AtomData(
                            atomInt1.exponent().get(i), atomInt2.exponent().get(j),
                            atomInt1.coord(), atomInt2.coord(), nucCoord.get(n),
                            atomInt1.cartAngular(), atomInt2.cartAngular()
                    );
                    resultado += -numAtomic.get(n)*normOrbital(atomInt1.exponent().get(i), atomInt1.cartAngular())
                            * normOrbital(atomInt2.exponent().get(j), atomInt2.cartAngular())
                            * atomInt1.coefficient().get(i) * atomInt2.coefficient().get(j)
                            * potInt(atData);
                }
            }
        }

        return resultado;
    }

    private static double potInt(AtomData atData) {

        // Distance squared between centers
        double dist2 = Math.pow(atData.coord1.x() - atData.coord2.x(), 2) +
                Math.pow(atData.coord1.y() - atData.coord2.y(), 2) +
                Math.pow(atData.coord1.z() - atData.coord2.z(), 2);

        double enerAB = Math.exp(-(atData.orbC1 * atData.orbC2 / (atData.orbC1 + atData.orbC2)) * dist2);

        return enerAB * 2.0 * Math.PI / (atData.orbC1 + atData.orbC2) * intChebyshev(atData);
    }

    public static double intChebyshev(
            AtomData atData) {
        final double PI = Math.PI;

        double error = 10.0;
        int n = 3;
        int j = 0;

        double varC0 = Math.cos(PI / 6.0);
        double varS0 = Math.sin(PI / 6.0);

        double varC1 = varS0;
        double varS1 = varC0;

        double q = (func1(abscissa(2, 1), atData)
                + func1(-abscissa(2, 1), atData)) * omega(2, 1);
        double p = func1(0.0, atData);
        double varCHP = p + q;

        while (error > Parameter.TOL_INTCHEVYCHEV && ((2 * n * (1 - j) + j * 4 * n / 3.0 - 1) <= Parameter.MAX_POINT_INTCHEVYCHEV)) {
            j = 1 - j;

            varC1 = j * varC1 + (1 - j) * varC0;
            varS1 = j * varS1 + (1 - j) * varS0;
            varC0 = j * varC0 + (1 - j) * Math.sqrt((1.0 + varC0) * 0.5);
            varS0 = j * varS0 + (1 - j) * (varS0 / (varC0 + varC0));

            double varC = varC0;
            double varS = varS0;

            // Inner loop
            for (int i = 1; i <= n - 1; i += 2) {
                double x = 1 + 2.0 / (3.0 * PI) * varS * varC * (3 + 2 * varS * varS) - ((double) i / (double) n);

                if (3 * Math.rint((i + j + j) / 3.0) > (i + j)) {
                    varCHP += (func1(-x, atData) +
                            func1(x, atData)) * Math.pow(varS, 4);
                }

                double temp = varS;
                varS = varS * varC1 + varC * varS1;
                varC = varC * varC1 - temp * varS1;
            }

            n = (1 + j) * n;
            p = p + (1 - j) * (varCHP - q);

            // Error estimate
            error = 16 * Math.abs((1 - j) * (q - 3 * p / 2.0) + j * (varCHP - 2 * q)) / (3.0 * n);
            q = (1 - j) * q + j * varCHP;
        }

        varCHP = 16 * q / (3 * n);
        return varCHP;
    }

    private static double func1(double x, AtomData atData) {
        double[] cons1 = new double[3];
        cons1[0] = ((atData.orbC1 * atData.coord1.x() + atData.orbC2 * atData.coord2.x()) / (atData.orbC1 + atData.orbC2)) - atData.coordNuc.x();
        cons1[1] = ((atData.orbC1 * atData.coord1.y() + atData.orbC2 * atData.coord2.y()) / (atData.orbC1 + atData.orbC2)) - atData.coordNuc.y();
        cons1[2] = ((atData.orbC1 * atData.coord1.z() + atData.orbC2 * atData.coord2.z()) / (atData.orbC1 + atData.orbC2)) - atData.coordNuc.z();

        AtomCoordData atXData = new AtomCoordData(
                atData.orbC1, atData.orbC2,
                atData.coord1.x(), atData.coord2.x(), atData.coordNuc.x()
        );

        AtomCoordData atYData = new AtomCoordData(
                atData.orbC1, atData.orbC2,
                atData.coord1.y(), atData.coord2.y(), atData.coordNuc.y()
        );

        AtomCoordData atZData = new AtomCoordData(
                atData.orbC1, atData.orbC2,
                atData.coord1.z(), atData.coord2.z(), atData.coordNuc.z()
        );

        return 0.5*Math.exp(-((atData.orbC1 + atData.orbC2) * Math.pow((x + 1) / 2,2) * dotProduct(cons1, cons1)))
                * potRecurr(atXData, atData.angCoord1.x(), atData.angCoord2.x(), x)
                * potRecurr(atYData, atData.angCoord1.y(), atData.angCoord2.y(), x)
                * potRecurr(atZData, atData.angCoord1.z(), atData.angCoord2.z(), x);
    }

    private static double potRecurr(AtomCoordData atCData, int angCoord1, int angCoord2, double x) {

        final double orbSum = atCData.orbC1 + atCData.orbC2;
        final double weightedCenter = (atCData.orbC1 * atCData.coord1 + atCData.orbC2 * atCData.coord2) / orbSum;
        final double t = (x + 1) / 2;

        // Case 1: both angular momenta are zero
        if (angCoord1 == 0 && angCoord2 == 0) {
            return 1.0;
        } else if (angCoord1 == 1 && angCoord2 == 0) {
            return -(atCData.coord1 - weightedCenter +
                    Math.pow(t,2)*(weightedCenter - atCData.nucCoord));
        } else if (angCoord1 > 1 && angCoord2 == 0) {
            return -(atCData.coord1 - weightedCenter +
                    Math.pow(t,2)*(weightedCenter - atCData.nucCoord))
                    * potRecurr(atCData, angCoord1 - 1, angCoord2, x)
                    +(angCoord1 - 1)/(2.0*(atCData.orbC1 + atCData.orbC2))*(1 - Math.pow(t, 2))
                    * potRecurr(atCData, angCoord1 -2, angCoord2, x);
        } else {
            return potRecurr(atCData, angCoord1 + 1, angCoord2 - 1, x)
                    + (atCData.coord1 - atCData.coord2)
                    * potRecurr(atCData, angCoord1, angCoord2 - 1, x);
        }

    }

    public static double normOrbital(double orbC, Atom.AngCoordinates vecA) {

        double pi = Math.PI;
        double b = Math.pow((2.0 * orbC / pi), 0.75);
        b *= Math.pow((4 * orbC), (vecA.x() + vecA.y() + vecA.z()) / 2.0);

        b /= Math.sqrt(
                doubleFactorial(2 * vecA.x() - 1)
                        * doubleFactorial(2 * vecA.y() - 1)
                        * doubleFactorial(2 * vecA.z() - 1)
        );

        return b;
    }

    public static double doubleFactorial(int n) {
        if (n <= 0) return 1.0;
        double result = 1.0;
        for (int i = n; i > 0; i -= 2) {
            result *= i;
        }
        return result;
    }

    public static double abscissa(int n, int i) {
        double theta = i * Math.PI / (n + 1.0);

        double term1 = (n + 1.0 - 2.0 * i) / (n + 1.0);
        double term2 = (2.0 / Math.PI) * (1.0 + (2.0 / 3.0) * Math.pow(Math.sin(theta), 2))
                * Math.sin(theta) * Math.cos(theta);
        return term1 + term2;
    }

    public static double omega(int n, int i) {
        int divInt = 16 / (3*(n + 1));
        return divInt * Math.pow(Math.sin((i * Math.PI) / (n + 1)), 4);
    }

    /**
     * Computes the dot product of two vectors of the same length.
     *
     * @param vecA first vector
     * @param vecB second vector
     * @return dot product vecA · vecB
     * @throws IllegalArgumentException if vectors have different lengths
     */
    public static double dotProduct(double[] vecA, double[] vecB) {
        if (vecA.length != vecB.length) {
            throw new IllegalArgumentException("Vectors must have the same length");
        }

        double sum = 0.0;
        for (int i = 0; i < vecA.length; i++) {
            sum += vecA[i] * vecB[i];
        }
        return sum;
    }

}
