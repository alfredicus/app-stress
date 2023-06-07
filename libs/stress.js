(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory(require("@youwol/math"));
	else if(typeof define === 'function' && define.amd)
		define("@alfredo-taboada/stress", ["@youwol/math"], factory);
	else if(typeof exports === 'object')
		exports["@alfredo-taboada/stress"] = factory(require("@youwol/math"));
	else
		root["@alfredo-taboada/stress"] = factory(root["@youwol/math"]);
})((typeof self !== 'undefined' ? self : this), (__WEBPACK_EXTERNAL_MODULE__youwol_math__) => {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "./lib/InverseMethod.ts":
/*!******************************!*\
  !*** ./lib/InverseMethod.ts ***!
  \******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "InverseMethod": () => (/* binding */ InverseMethod),
/* harmony export */   "cloneMisfitCriteriunSolution": () => (/* binding */ cloneMisfitCriteriunSolution),
/* harmony export */   "createDefaultSolution": () => (/* binding */ createDefaultSolution)
/* harmony export */ });
/* harmony import */ var _search_GridSearch__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./search/GridSearch */ "./lib/search/GridSearch.ts");
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./types/math */ "./lib/types/math.ts");


/**
 * @category Inversion
 */
function cloneMisfitCriteriunSolution(misfitCriteriunSolution) {
    return {
        // criterion: misfitCriteriunSolution.criterion,
        misfit: misfitCriteriunSolution.misfit,
        rotationMatrixW: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(misfitCriteriunSolution.rotationMatrixW),
        rotationMatrixD: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(misfitCriteriunSolution.rotationMatrixD),
        stressRatio: misfitCriteriunSolution.stressRatio,
        stressTensorSolution: misfitCriteriunSolution.stressTensorSolution
    };
}
/**
 * @category Inversion
 */
function createDefaultSolution() {
    return {
        misfit: Number.POSITIVE_INFINITY,
        rotationMatrixW: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)(),
        rotationMatrixD: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)(),
        stressRatio: 0,
        stressTensorSolution: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)()
    };
}
/**
 * @category Inversion
 */
class InverseMethod {
    constructor() {
        this.misfitCriteriunSolution = {
            misfit: Number.POSITIVE_INFINITY,
            rotationMatrixW: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)(),
            rotationMatrixD: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)(),
            stressRatio: 0,
            stressTensorSolution: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)()
        };
        this.searchMethod = new _search_GridSearch__WEBPACK_IMPORTED_MODULE_0__.GridSearch();
        this.data_ = [];
    }
    get data() {
        return this.data_;
    }
    setSearchMethod(search) {
        this.searchMethod = search;
    }
    addData(data) {
        if (Array.isArray(data)) {
            data.forEach(d => this.data_.push(d));
        }
        else {
            this.data_.push(data);
        }
    }
    run(reset = true) {
        if (this.data_.length === 0) {
            throw new Error('No data provided');
        }
        if (reset) {
            this.misfitCriteriunSolution.misfit = Number.POSITIVE_INFINITY;
        }
        return this.searchMethod.run(this.data_, this.misfitCriteriunSolution);
    }
    cost({ displ, strain, stress }) {
        if (this.data_.length === 0) {
            throw new Error('No data provided');
        }
        return this.data_.reduce((cumul, data) => cumul + data.cost({ displ, strain, stress }), 0) / this.data_.length;
    }
}


/***/ }),

/***/ "./lib/analysis/Curve3D.ts":
/*!*********************************!*\
  !*** ./lib/analysis/Curve3D.ts ***!
  \*********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "Curve3D": () => (/* binding */ Curve3D)
/* harmony export */ });
class Curve3D {
    constructor() {
        this.points = [];
    }
    addPoint(x, y, z = 0) {
        if (typeof x === 'number') {
            this.points.push(x, y, z);
        }
        else {
            this.points.push(x[0], x[1], x[2]);
        }
    }
    clear() {
        this.points = [];
    }
    get buffer() {
        let s = `GOCAD PLine 1.0
        HEADER {
            name: curve3d
        }
        PROPERTIES rake
        `;
        const nbPoints = this.points.length / 3;
        let index = 0;
        for (let i = 0; i < this.points.length; i += 3) {
            const attr = 0;
            s += 'PVRTX ' + index + ' ' + this.points[i] + ' ' + this.points[i + 1] + ' ' + this.points[i + 2] + ' ' + attr + '\n';
            index++;
        }
        for (let i = 0; i < nbPoints - 1; ++i) {
            s += 'SEG ' + i + ' ' + (i + 1) + '\n';
        }
        s += 'END';
        return s;
    }
}


/***/ }),

/***/ "./lib/analysis/EquipotentialCurve.ts":
/*!********************************************!*\
  !*** ./lib/analysis/EquipotentialCurve.ts ***!
  \********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "EquipotentialCurve": () => (/* binding */ EquipotentialCurve)
/* harmony export */ });
/* harmony import */ var _types_mechanics__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types/mechanics */ "./lib/types/mechanics.ts");
/* harmony import */ var _Curve3D__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Curve3D */ "./lib/analysis/Curve3D.ts");


class EquipotentialCurve {
    constructor(lambda, radius = 1) {
        this.lambda = lambda;
        this.r = radius;
        this.exp_sin = (this.lambda[0] - this.lambda[2]) / (this.lambda[1] - this.lambda[0]);
        this.exp_cos = (this.lambda[2] - this.lambda[1]) / (this.lambda[1] - this.lambda[0]);
    }
    generate(theta, phi) {
        const lineBuilder = new _Curve3D__WEBPACK_IMPORTED_MODULE_1__.Curve3D();
        theta = theta * Math.PI / 180;
        phi = phi * Math.PI / 180;
        // NORMAL STRESS EQUIPOTENTIALS 
        // The principal stress values are defined according to the rock mechanics sign convention (positive values for compressive stresses)
        const sigma_1 = -this.lambda[0]; // Principal stress in X direction
        const sigma_2 = -this.lambda[2]; // Principal stress in Z direction
        const sigma_3 = -this.lambda[1]; // Principal stress in Y direction
        // Center and radius of Mohr circle 3 - 1
        const sigma_3_1_mean = (sigma_1 + sigma_3) / 2;
        const rad_3_1 = (sigma_1 - sigma_3) / 2;
        // Center and radius of Mohr circle 3 - 2
        const sigma_3_2_mean = (sigma_2 + sigma_3) / 2;
        const rad_3_2 = (sigma_2 - sigma_3) / 2;
        // Center and radius of Mohr circle 2 - 1
        const sigma_2_1_mean = (sigma_1 + sigma_2) / 2;
        const rad_2_1 = (sigma_1 - sigma_2) / 2;
        // The integral lines derive from a scalar function defined by the normal stress component sigma_n
        // sigma_n is calculated for a specific plane tangent to the sphere whose normal vector is defined by (phi_1, theta_1)
        // sigma_n = (this.lambda[0] * cos(phi)**2 +this.lambda[1] * sin(phi)**2 ) * sin(theta)**2 + this.lambda[2] * ( 1 - sin(theta)**2 )
        let sigma_n = (this.lambda[0] * Math.cos(phi) ** 2 + this.lambda[1] * Math.sin(phi) ** 2) * Math.sin(theta) ** 2 + this.lambda[2] * Math.cos(theta) ** 2;
        sigma_n *= -1;
        // Plot equipotential corresponding to the normal force that passes through the fixed point:
        //      
        // tau_n0 = shear stress for point pO located in the Mohr circle 3-1 at normal stress sigma_n
        const tau_0 = Math.sqrt(rad_3_1 ** 2 - (sigma_n - sigma_3_1_mean) ** 2);
        // alfa_n0 = angle (in radians) between line segment (sigma_3,tau_n0) in circle 3-1 and the horizontal axis (sigma_3, sigma_1) :
        const alfa_0 = Math.atan(tau_0 / (sigma_n - sigma_3));
        /* if ( sigma_2 == sigma_3) {
            // Particular Case 1: revolution stress tensor around sigma_1
            // Equipotential curve is a circle sweeping at an angle alfa_n0 around sigma_1
            arcCircle(r: this.r, 'sigma_1', alfa_0 )
        }
        else if ( sigma_2 == sigma_1) {
            // Particular Case 2: revolution stress tensor around sigma_3
            // Equipotential curve is a circle sweeping at an angle PI/2 - alfa_n0 around sigma_3
            arcCircle(r: this.r, 'sigma_3', Math.PI/2 - alfa_0 )
        } */
        if (sigma_n > sigma_3 && sigma_n < sigma_2) {
            // Case 1: the equipotential line lies between circle 3 - 1 and circle 3 - 2:
            // tau_1 = shear stress for point p1 located in the Mohr circle 3-2 at normal stress sigma_n
            const tau_1 = Math.sqrt(rad_3_2 ** 2 - (sigma_n - sigma_3_2_mean) ** 2);
            // alfa_1 = angle (in radians) between line segment (sigma_3,tau_n0) in circle 3-2 and the horizontal axis (sigma_3, sigma_2) :
            const alfa_1 = Math.atan(tau_1 / (sigma_n - sigma_3));
            // Plot curve corresponding to the line segment between points: (sigma_n, tau_0) and (sigma_n, tau_1)
            const first = {
                circle: '3_1',
                p: [sigma_n, tau_0],
                angle: alfa_0
            };
            const second = {
                circle: '3_2',
                p: [sigma_n, tau_1],
                angle: alfa_1
            };
            return (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_0__.mohrCircleLine)({ r: this.r, first, second, sigma_1, sigma_2, sigma_3 });
        }
        else if (sigma_n > sigma_2 && sigma_n < sigma_1) {
            // Case 2: the equipotential line lies between circle 3 - 1 and circle 2 - 1:
            // tau_1 = shear stress for point pO located in the Mohr circle 3-1 at normal stress sigma_n
            const tau_1 = Math.sqrt(rad_2_1 ** 2 - (sigma_n - sigma_2_1_mean) ** 2);
            // alfa_1 = angle (in radians) between line segment (sigma_2,tau_n0) in circle 2-1 and the horizontal axis (sigma_2, sigma_1) :
            const alfa_1 = Math.atan(tau_1 / (sigma_n - sigma_2));
            // Plot curve corresponding to the line segment between points: (sigma_n, tau_n0) and (sigma_n, tau_1)
            const first = {
                circle: '3_1',
                p: [sigma_n, tau_0],
                angle: alfa_0
            };
            const second = {
                circle: '2_1',
                p: [sigma_n, tau_1],
                angle: alfa_1
            };
            return (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_0__.mohrCircleLine)({ r: this.r, first, second, sigma_1, sigma_2, sigma_3 });
        }
        else if (sigma_n == sigma_2) {
            // Case 3: the equipotential line lies between circle 3 - 1 and sigma_2:
            // tau_1 = shear stress for point pO
            const tau_1 = 0;
            // alfa_1 = angle (in radians) between line segment (sigma_3,tau_n0) in circle 3-2 and the horizontal axis (sigma_3, sigma_2) :
            const alfa_1 = 0;
            // Plot curve corresponding to the line segment between points: (sigma_n, tau_n0) and (sigma_n, tau_1)
            const first = {
                circle: '3_1',
                p: [sigma_n, tau_0],
                angle: alfa_0
            };
            const second = {
                circle: '3_2',
                p: [sigma_n, tau_1],
                angle: alfa_1
            };
            return (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_0__.mohrCircleLine)({ r: this.r, first, second, sigma_1, sigma_2, sigma_3 });
        }
        return '';
    }
}


/***/ }),

/***/ "./lib/analysis/IntegralCurve.ts":
/*!***************************************!*\
  !*** ./lib/analysis/IntegralCurve.ts ***!
  \***************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "IntegralCurve": () => (/* binding */ IntegralCurve)
/* harmony export */ });
/* harmony import */ var _Curve3D__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Curve3D */ "./lib/analysis/Curve3D.ts");

/**
 * @brief Calculate the integral curves on the sphere surface that are parallel to the shear stress (streamlines in fluid mechanics)
 * @param {r: number, exp_cos: number, exp_sin: number}
 * @example
 * ```ts
 *
 * ```
 */
class IntegralCurve {
    constructor(lambda, radius = 1) {
        this.lambda = lambda;
        this.r = radius;
        this.exp_sin = (this.lambda[0] - this.lambda[2]) / (this.lambda[1] - this.lambda[0]);
        this.exp_cos = (this.lambda[2] - this.lambda[1]) / (this.lambda[1] - this.lambda[0]);
    }
    /**
     * Define a tangent plane by fixing angles phi and theta in spherical coordinates
     * @param theta In degrees
     * @param phi In degrees
     */
    generate(theta, phi) {
        const lineBuilder = new _Curve3D__WEBPACK_IMPORTED_MODULE_0__.Curve3D();
        let phi_1 = Math.PI * phi / 180;
        let theta_1 = Math.PI * theta / 180;
        // Determine integral curve that passes by this point of the sphere surface, by calculating k1
        // k1 is defined as a function of phi and theta for a specific symmetrical tensor
        let k1 = Math.tan(theta_1) / (Math.sin(phi_1) ** this.exp_sin * Math.cos(phi_1) ** this.exp_cos);
        // Plot the integral curve that passes by this specific point
        for (let i = 0; i <= 180; ++i) {
            const phi = Math.PI * i / 360;
            const theta = Math.atan(k1 * Math.sin(phi) ** this.exp_sin * Math.cos(phi) ** this.exp_cos);
            const x = this.r * Math.sin(theta) * Math.cos(phi);
            const y = this.r * Math.sin(theta) * Math.sin(phi);
            const z = this.r * Math.cos(theta);
            lineBuilder.addPoint(x, y, z);
        }
        return lineBuilder.buffer;
    }
}


/***/ }),

/***/ "./lib/analysis/MohrCoulombCurves.ts":
/*!*******************************************!*\
  !*** ./lib/analysis/MohrCoulombCurves.ts ***!
  \*******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "MohrCoulombCurve": () => (/* binding */ MohrCoulombCurve)
/* harmony export */ });
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");
/* harmony import */ var _types_mechanics__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types/mechanics */ "./lib/types/mechanics.ts");
/* harmony import */ var _Curve3D__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Curve3D */ "./lib/analysis/Curve3D.ts");



class MohrCoulombCurve {
    constructor(lambda, radius = 1) {
        this.lambda = lambda;
        this.r = radius;
    }
    maxFrictionAngle(cohesion) {
        const sigma_1 = -this.lambda[0]; // Principal stress in X direction
        const sigma_3 = -this.lambda[1]; // Principal stress in Y direction
        // Center and radius of Mohr circle 3 - 1
        const sigma_3_1_mean = (sigma_1 + sigma_3) / 2;
        const rad_3_1 = (sigma_1 - sigma_3) / 2;
        // The friction angle phi_a varies in the range [0, phi_3_1], 
        //      where phi_3_1 is the angle of the line tangent to Mohr circle between sigma_3 - sigma_1 that passes by c
        let a_3_1 = sigma_3_1_mean ** 2 + cohesion ** 2;
        let b_3_1 = -2 * rad_3_1 * sigma_3_1_mean;
        let c_3_1 = rad_3_1 ** 2 - cohesion ** 2;
        // Angle max is defined in degrees
        const max = Math.asin((-b_3_1 - Math.sqrt(b_3_1 ** 2 - 4 * a_3_1 * c_3_1)) / (2 * a_3_1)) * 180 / Math.PI;
        return max;
    }
    maxCohesion(frictionAngle) {
        const sigma_1 = -this.lambda[0]; // Principal stress in X direction
        const sigma_3 = -this.lambda[1]; // Principal stress in Y direction
        // Center and radius of Mohr circle 3 - 1
        const sigma_3_1_mean = (sigma_1 + sigma_3) / 2;
        const rad_3_1 = (sigma_1 - sigma_3) / 2;
        // The friction angle phi_a varies in the range [0, phi_3_1], 
        //      where phi_3_1 is the angle of the line tangent to Mohr circle between sigma_3 - sigma_1 that passes by c
        let frictionAngleRad = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(frictionAngle);
        const cohesion = (rad_3_1 - sigma_3_1_mean * Math.sin(frictionAngleRad)) / Math.cos(frictionAngleRad);
        return cohesion;
    }
    /**
     * Define a tangent plane by fixing angles phi and theta in spherical coordinates
     * @param frictionAngle Friction angle in degrees
     * @param cohesion In [0, 0.5]
     */
    generate(frictionAngle, cohesion) {
        const lineBuilder = new _Curve3D__WEBPACK_IMPORTED_MODULE_2__.Curve3D();
        // The principal stress values are defined according to the rock mechanics sign convention (positive values for compressive stresses)
        const sigma_1 = -this.lambda[0]; // Principal stress in X direction
        const sigma_2 = -this.lambda[2]; // Principal stress in Z direction
        const sigma_3 = -this.lambda[1]; // Principal stress in Y direction
        // Center and radius of Mohr circle 3 - 1
        const sigma_3_1_mean = (sigma_1 + sigma_3) / 2;
        const rad_3_1 = (sigma_1 - sigma_3) / 2;
        // Center and radius of Mohr circle 3 - 2
        const sigma_3_2_mean = (sigma_2 + sigma_3) / 2;
        const rad_3_2 = (sigma_2 - sigma_3) / 2;
        // Center and radius of Mohr circle 2 - 1
        const sigma_2_1_mean = (sigma_1 + sigma_2) / 2;
        const rad_2_1 = (sigma_1 - sigma_2) / 2;
        if (cohesion < 0 || cohesion > 0.5) {
            throw new Error('Cohesion must be in [0, 0.5]');
        }
        const c = cohesion;
        const phi_a = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(frictionAngle);
        // MOHR-COULOMB FRICTION
        // The frictional strength law is defined by a Mohr-Coulomb line involving friction and cohesion:
        // These two parameters must be fixed by the user in the menu within a predefined range
        // The cohesion c is defined between [0,0.5] for a normalized stress tensor (sigma_1=1, sigma_2= R, sigma_3=0)
        // The friction angle phi_a varies in the range [0, phi_3_1], 
        //      where phi_3_1 is the angle of the line tangent to Mohr circle between sigma_3 - sigma_1 that passes by c        
        let a_3_1 = sigma_3_1_mean ** 2 + c ** 2;
        let b_3_1 = -2 * rad_3_1 * sigma_3_1_mean;
        let c_3_1 = rad_3_1 ** 2 - c ** 2;
        // Angles are defined in radians
        const phi_3_1 = Math.asin((-b_3_1 - Math.sqrt(b_3_1 ** 2 - 4 * a_3_1 * c_3_1)) / (2 * a_3_1));
        // We calculate two other threshold angles for lines tangent to Mohr circles between sigma_3 - sigma_2, and sigma_2 - sigma_1
        // Angle phi_3_2 for Mohr circle between sigma_2 - sigma_3 is calculated with a similar equation:
        //      phi_3_2 is calculated only for positive values
        let phi_3_2_bool = false; // Before verification, we suppose that the friction line does not intersect Mohr Circle S3-S2
        let phi_3_2 = 0;
        if (c < rad_3_2) {
            phi_3_2_bool = true;
            // phi_3_2 > 0, thus the friction line may intersect Mohr circle between sigma_3 - sigma_2
            // phi_3_2 can be calculated from a trigonometric expression
            let a_3_2 = sigma_3_2_mean ** 2 + c ** 2;
            let b_3_2 = -2 * rad_3_2 * sigma_3_2_mean;
            let c_3_2 = rad_3_2 ** 2 - c ** 2;
            phi_3_2 = Math.asin((-b_3_2 - Math.sqrt(b_3_2 ** 2 - 4 * a_3_2 * c_3_2)) / (2 * a_3_2));
        }
        // Angle phi_2_1 for Mohr circle between sigma_1 - sigma_2 is calculated from the general equation:
        //      tan(pi / 4 + phi_2_1 / 2) = ( - c + sqrt( c^2 + sigma_2 * sigma_1 ) ) / sigma_2
        //      phi_2_1 is calculated only for positive values
        let phi_2_1_bool = false; // Before verification, we suppose that the friction line does not intersect Mohr Circle S2-S1
        let phi_2_1 = 0;
        if (c < rad_2_1) {
            phi_2_1_bool = true;
            // phi_2_1 > 0, thus the friction line may intersect Mohr circle between sigma_3 - sigma_2
            // phi_2_1 can be calculated from a trigonometric expression
            let a_2_1 = sigma_2_1_mean ** 2 + c ** 2;
            let b_2_1 = -2 * rad_2_1 * sigma_2_1_mean;
            let c_2_1 = rad_2_1 ** 2 - c ** 2;
            phi_2_1 = Math.asin((-b_2_1 - Math.sqrt(b_2_1 ** 2 - 4 * a_2_1 * c_2_1)) / (2 * a_2_1));
        }
        // Let (c, phi_a) be the cohesion and friction angle chosen by the user. The friction coefficient mu is the slope of the line
        const mu = Math.tan(phi_a);
        const sigma_3_1 = [0, 0];
        const tau_3_1 = [0, 0];
        const alfa_3_1 = [0, 0];
        const sigma_3_2 = [0, 0];
        const tau_3_2 = [0, 0];
        const alfa_3_2 = [0, 0];
        const sigma_2_1 = [0, 0];
        const tau_2_1 = [0, 0];
        const alfa_2_1 = [0, 0];
        const struct = {
            sigma_3_1,
            tau_3_1,
            alfa_3_1,
            sigma_3_2,
            tau_3_2,
            alfa_3_2,
            sigma_2_1,
            tau_2_1,
            alfa_2_1,
        };
        // The Mohr-Coulomb line intersects the 3D Mohr circle in 1, 2 or 3 different line segments 
        // The line always intersects circle 1 - 3 in two points named: sigma_3_1[0], tau_3_1[0] and sigma_3_1[1], tau_3_1[1]
        //      such that sigma_3_1[0] <= sigma_3_1[1]
        // sigma_3_1 values are given by the roots of a quadratic equation a x^2 + b x + c = 0, with coeffcients a, b, anc c, as follows:
        const a_qe = 1 + mu ** 2;
        const b_qe = 2 * (c * mu - sigma_3_1_mean);
        const c_qe = c ** 2 + sigma_1 * sigma_3;
        //  Calculate the discriminant of the quadratic equation:
        let delta = Math.sqrt(b_qe ** 2 - 4 * a_qe * c_qe);
        //  Calculate intersection points
        sigma_3_1[0] = (-b_qe - delta) / (2 * a_qe);
        tau_3_1[0] = Math.sqrt(rad_3_1 ** 2 - (sigma_3_1[0] - sigma_3_1_mean) ** 2);
        sigma_3_1[1] = (-b_qe + delta) / (2 * a_qe);
        tau_3_1[1] = Math.sqrt(rad_3_1 ** 2 - (sigma_3_1[1] - sigma_3_1_mean) ** 2);
        // Calculate the angle (in radians) between segment (sigma_3,tau_3_1) in circle 3-1 and the horizontal axis (sigma_3, sigma_1) :
        // alfa_3_1 is the azimuthal angle in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
        alfa_3_1[0] = Math.atan(tau_3_1[0] / (sigma_3_1[0] - sigma_3));
        alfa_3_1[1] = Math.atan(tau_3_1[1] / (sigma_3_1[1] - sigma_3));
        // Define booleans indicating if the friction line intersects the smaller Mohr circles:
        let circle_2_1 = false;
        let circle_3_2 = false;
        if (phi_3_2_bool) {
            // phi_3_2 > 0, thus the friction line may intersect Mohr circle between sigma_3 - sigma_2
            if (phi_a < phi_3_2) {
                circle_3_2 = true;
                // The Mohr-Coulomb line intersects circle 2 - 3 in two points named: sigma_3_2[0], tau_3_2[0] and sigma_3_2[1], tau_3_2[1]
                //      such that sigma_3_2[0] <= sigma_3_2[1]
                // sigma_3_2 values are given by the roots of a quadratic equation a x^2 + b x + c = 0, with coeffcients a, b, anc c, as follows:
                const a_qe = 1 + mu ** 2;
                const b_qe = 2 * (c * mu - sigma_3_2_mean);
                const c_qe = c ** 2 + sigma_2 * sigma_3;
                //  Calculate the discriminant of the quadratic equation:
                let delta = Math.sqrt(b_qe ** 2 - 4 * a_qe * c_qe);
                //  Calculate intersection points
                sigma_3_2[0] = (-b_qe - delta) / (2 * a_qe);
                tau_3_2[0] = Math.sqrt(rad_3_2 ** 2 - (sigma_3_2[0] - sigma_3_2_mean) ** 2);
                sigma_3_2[1] = (-b_qe + delta) / (2 * a_qe);
                tau_3_2[1] = Math.sqrt(rad_3_2 ** 2 - (sigma_3_2[1] - sigma_3_2_mean) ** 2);
                // Calculate the angle (in radians) between segment (sigma_3,tau_3_2) in circle 3-2 and the horizontal axis (sigma_3, sigma_2) :
                // alfa_3_2 is the polar angle in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
                alfa_3_2[0] = Math.atan(tau_3_2[0] / (sigma_3_2[0] - sigma_3));
                alfa_3_2[1] = Math.atan(tau_3_2[1] / (sigma_3_2[1] - sigma_3));
            }
        }
        if (phi_2_1_bool) {
            // phi_2_1 > 0, thus the friction line may intersect Mohr circle between sigma_3 - sigma_2
            if (phi_a < phi_2_1) {
                circle_2_1 = true;
                // The Mohr-Coulomb line intersects circle 1 - 2 in two points named: sigma_2_1[0], tau_2_1[0] and sigma_2_1[1], tau_2_1[1]
                //      such that sigma_2_1[0] <= sigma_2_1[1]
                // sigma_2_1 values are given by the roots of a quadratic equation a x^2 + b x + c = 0, with coeffcients a, b, anc c, as follows:
                const a_qe = 1 + mu ** 2;
                const b_qe = 2 * (c * mu - sigma_2_1_mean);
                const c_qe = c ** 2 + sigma_1 * sigma_2;
                //  Calculate the discriminant of the quadratic equation:
                let delta = Math.sqrt(b_qe ** 2 - 4 * a_qe * c_qe);
                //  Calculate intersection points
                sigma_2_1[0] = (-b_qe - delta) / (2 * a_qe);
                tau_2_1[0] = Math.sqrt(rad_2_1 ** 2 - (sigma_2_1[0] - sigma_2_1_mean) ** 2);
                sigma_2_1[1] = (-b_qe + delta) / (2 * a_qe);
                tau_2_1[1] = Math.sqrt(rad_2_1 ** 2 - (sigma_2_1[1] - sigma_2_1_mean) ** 2);
                // Calculate the angle (in radians) between segment (sigma_2,tau_2_1) in circle 2-1 and the horizontal axis (sigma_2, sigma_1) :
                // alfa_2_1 is the latitude in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
                alfa_2_1[0] = Math.atan(tau_2_1[0] / (sigma_2_1[0] - sigma_2));
                alfa_2_1[1] = Math.atan(tau_2_1[1] / (sigma_2_1[1] - sigma_2));
            }
        }
        // We calculate the corresponding curves in the sphere:
        if (!circle_2_1 && !circle_3_2) {
            // Case 1: the Mohr-Coulomb line only intersects circle 3 - 1
            // Plot curve corresponding to the line segment between points: 
            //      sigma_3_1[0], tau_3_1[0] and sigma_3_1[1], tau_3_1[1]
            return (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(0, '3_1', struct),
                second: this.getPoint(1, '3_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            // return mohrCircleLine( {r: this.r, first, second, sigma_1, sigma_2, sigma_3} )
        }
        else if (!circle_2_1 && circle_3_2) {
            // Case 2: the Mohr-Coulomb line intersects circle 3 - 1 and circle 3 - 2
            // Plot curves corresponding to the line segment between points: 
            //      sigma_3_1[0], tau_3_1[0] and sigma_3_2[0], tau_3_2[0];
            let buffer = (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(0, '3_1', struct),
                second: this.getPoint(0, '3_2', struct),
                sigma_1, sigma_2, sigma_3
            });
            buffer += '\n';
            //      sigma_3_2[1], tau_3_2[1] and sigma_3_1[1], tau_3_1[1];
            buffer += (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(1, '3_2', struct),
                second: this.getPoint(1, '3_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            return buffer;
        }
        else if (circle_2_1 && !circle_3_2) {
            // Case 3: the Mohr-Coulomb line intersects circle 3 - 1 and circle 2 - 1
            // Plot curves corresponding to the line segment between points: 
            //      sigma_3_1[0], tau_3_1[0] and sigma_2_1[0], tau_2_1[0];
            let buffer = (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(0, '3_1', struct),
                second: this.getPoint(0, '2_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            buffer += '\n';
            //      sigma_2_1[1], tau_2_1[1] and sigma_3_1[1], tau_3_1[1];
            buffer += (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(1, '2_1', struct),
                second: this.getPoint(1, '3_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            return buffer;
        }
        else {
            // Case 4: the Mohr-Coulomb line intersects circle 3 - 1, circle 3 - 2, and circle 2 - 1
            // Plot curves corresponding to the line segment between points: 
            //      sigma_3_1[0], tau_3_1[0] and sigma_3_2[0], tau_3_2[0];
            let buffer = (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(0, '3_1', struct),
                second: this.getPoint(0, '3_2', struct),
                sigma_1, sigma_2, sigma_3
            });
            buffer += '\n';
            //      sigma_3_2[1], tau_3_2[1] and sigma_2_1[0], tau_2_1[0];
            buffer += (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(1, '3_2', struct),
                second: this.getPoint(0, '2_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            buffer += '\n';
            //      sigma_2_1[1], tau_2_1[1] and sigma_3_1[1], tau_3_1[1];
            buffer += (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_1__.mohrCircleLine)({
                r: this.r,
                first: this.getPoint(1, '2_1', struct),
                second: this.getPoint(1, '3_1', struct),
                sigma_1, sigma_2, sigma_3
            });
            return buffer;
        }
    }
    getPoint(index, name, struct) {
        if (name === '3_1') {
            return {
                circle: name,
                p: [struct.sigma_3_1[index], struct.tau_3_1[index]],
                angle: struct.alfa_3_1[index]
            };
        }
        else if (name === '3_2') {
            return {
                circle: name,
                p: [struct.sigma_3_2[index], struct.tau_3_2[index]],
                angle: struct.alfa_3_2[index]
            };
        }
        else if (name === '2_1') {
            return {
                circle: name,
                p: [struct.sigma_2_1[index], struct.tau_2_1[index]],
                angle: struct.alfa_2_1[index]
            };
        }
        else {
            throw new Error(`name ${name} is unknown. Should be 3_1, 3_2 or 2_1`);
        }
    }
}


/***/ }),

/***/ "./lib/analysis/index.ts":
/*!*******************************!*\
  !*** ./lib/analysis/index.ts ***!
  \*******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "Curve3D": () => (/* reexport safe */ _Curve3D__WEBPACK_IMPORTED_MODULE_0__.Curve3D),
/* harmony export */   "EquipotentialCurve": () => (/* reexport safe */ _EquipotentialCurve__WEBPACK_IMPORTED_MODULE_1__.EquipotentialCurve),
/* harmony export */   "IntegralCurve": () => (/* reexport safe */ _IntegralCurve__WEBPACK_IMPORTED_MODULE_2__.IntegralCurve),
/* harmony export */   "MohrCoulombCurve": () => (/* reexport safe */ _MohrCoulombCurves__WEBPACK_IMPORTED_MODULE_3__.MohrCoulombCurve)
/* harmony export */ });
/* harmony import */ var _Curve3D__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Curve3D */ "./lib/analysis/Curve3D.ts");
/* harmony import */ var _EquipotentialCurve__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./EquipotentialCurve */ "./lib/analysis/EquipotentialCurve.ts");
/* harmony import */ var _IntegralCurve__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./IntegralCurve */ "./lib/analysis/IntegralCurve.ts");
/* harmony import */ var _MohrCoulombCurves__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./MohrCoulombCurves */ "./lib/analysis/MohrCoulombCurves.ts");






/***/ }),

/***/ "./lib/data/CompactionShearBands_Kin.ts":
/*!**********************************************!*\
  !*** ./lib/data/CompactionShearBands_Kin.ts ***!
  \**********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "CompactionShearBandsKin": () => (/* binding */ CompactionShearBandsKin)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _Data__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Data */ "./lib/data/Data.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./types */ "./lib/data/types.ts");
/* harmony import */ var _utils_CompactionShearBands__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/CompactionShearBands */ "./lib/utils/CompactionShearBands.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils */ "./lib/utils/index.ts");
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");
// import {
//     Matrix3x3, Point3D, scalarProductUnitVectors,
//     SphericalCoords, Vector3, add_Vectors,
//     constant_x_Vector, deg2rad, spherical2unitVectorCartesian,
//     faultStressComponents, normalizeVector, unitVectorCartesian2Spherical,
//     normalizedCrossProduct,
//     multiplyTensors
// }
//     from "../types"
// import { Data, DataParameters } from "./Data"
// import { StriatedPlaneProblemType } from "./StriatedPlane_Kin"
// import { FractureStrategy } from "./types"
// import { Fault, Direction, SensOfMovement } from "../utils/Fault"
// // import { Fracture, FractureParams, FractureStrategy } from "./Fracture"


// import { faultStressComponents } from "../types/mechanics"
// import { Fault, Direction, SensOfMovement } from "../utils/Fault"




/**
 Compaction-Shear Bands:
 A pair of compaction-shear bands is defined by two planes whose plane of movement is perpendicular to the intersection line between the planes.
 The plane of movement is defined by the two normal vectors to the compaction-shear bands.
 We make the following hypotheses concerning principal stress orientations:
    The compressional axis Sigma 1 is located in the plane of movement and bisects the obtuse angle (>= 90°) between compaction-shear bands
    The extensional axis Sigma 3 is located in the plane of movement and bisects the acute angle (<= 90°) between compaction-shear bands

 Compaction-shear bands are defined in the input file in two consecutive lines.
 Each line specifies all the available data for each compaction-shear band.

 Several data sets defining two compaction-shear bands are considered:
 1) Case 1: The geometry and kinematics of the compaction-shear bands are defined, yet the rake angles are NOT defined:
    The orientation of the principal axes are calculated from the geometry of the compaction-shear bands.
    The intermediate axis Sigma 2 is parallel to the intersection line between the compaction-shear bands;
    The intermediate axis Sigma 2 is perpendicular to the plane of mouvement;
    a) Each plane is defined by a set of 3 parameters, as follows:
        Fault strike: clockwise angle measured from the North direction [0, 360)
        Fault dip: [0, 90]
        Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
    b) The rake angles defining the slip vectors are NOT defined
    c) The sense of mouvement is indicated for each compaction-shear band :
        For verification purposes, it is recommended to indicate both the dip-slip and strike-slip compoenents, when possible.
          Dip-slip component:
              N = Normal fault,
              I = Inverse fault or thrust
          Strike-slip componenet:
              RL = Right-Lateral fault
              LL = Left-Lateral fault
        Sense of mouvement: N, I, RL, LL, N-RL, N-LL, I-RL, I-LL
2) Case 2: The geometry, the striation (i.e., the rake), and the kinematics of one or both compaction-shear bands are defined.
    This case is not considered as the striations define redundant information for constraining the stress tensor orientation.
    Compaction-shear bands with striations can be considered separately as neoformed striated planes for which the angle <Sigma 1,n> is known


3) Particular case:
    If the striation is known for one of the conjugate faults, then the other plane
    can be defined by one axis that is contained in the plane.
    However, this case would require a special format in the input file, which is inconvenient...

    @category Data
 */
class CompactionShearBandsKin extends _Data__WEBPACK_IMPORTED_MODULE_1__.Data {
    constructor() {
        super(...arguments);
        this.nPlane = [];
        this.nStriation1 = undefined;
        this.nStriation2 = undefined;
        this.pos = undefined;
        this.problemType = _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC;
        this.strategy = _types__WEBPACK_IMPORTED_MODULE_2__.FractureStrategy.ANGLE;
        this.oriented = true;
        this.EPS = 1e-7;
        this.plane = [false, false];
        // Principal directions of the data
        this.nSigma1_Sc = undefined;
        this.nSigma2_Sc = undefined;
        this.nSigma3_Sc = undefined;
        this.Crot = undefined;
        this.csb = [undefined, undefined];
        this.params = [undefined, undefined];
        this.striation = [true, true];
    }
    nbLinkedData() {
        return 2;
    }
    initialize(params) {
        let nPlane2Neg;
        // TODO
        // CHeck that params.length === 2
        for (let i = 0; i < 2; ++i) {
            if (Number.isNaN(params[i].strike)) {
                throw new Error('Missing strike angle for compaction-shear band ' + params[i].noPlane);
            }
            if (Number.isNaN(params[i].dip)) {
                throw new Error('Missing dip angle for compaction-shear band ' + params[i].noPlane);
            }
            if (params[i].dip < 90 && params[i].dipDirection === undefined) {
                throw new Error('Missing dip direction for compaction-shear band  ' + params[i].noPlane);
            }
            this.params[i] = {
                noPlane: params[i].noPlane,
                strike: params[i].strike,
                dipDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[i].dipDirection),
                dip: params[i].dip,
                sensOfMovement: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getSensOfMovementFromString)(params[i].typeOfMovement),
                rake: params[i].rake,
                strikeDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[i].strikeDirection)
            };
            this.csb[i] = _utils_CompactionShearBands__WEBPACK_IMPORTED_MODULE_3__.CompactionShearBands.create(this.params[i]);
            // compaction-shear band i is defined: (strike1, dip1, dipDirection1)
            this.plane[i] = true;
            // Calculate the unit vector normal to plane i: nPlane
            this.nPlane[i] = this.csb[i].nPlane;
        }
        if (this.nPlane[0] === this.nPlane[1] || this.nPlane[0] === (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -1, V: this.nPlane[1] })) {
            throw new Error('The two compaction-shear bands ' + this.params[0].noPlane + ' and ' + this.params[1].noPlane + ' are identical');
        }
        /** this.nStriation = nStriation
        // this.nPerpStriation = nPerpStriation

        // Check orthogonality
        // const sp = scalarProductUnitVectors({U: nPlane, V: nStriation})
        //*if (Math.abs(sp) >this.EPS) {
            throw new Error(`striation is not on the fault plane. Dot product gives ${sp}`)
        } **/
        return true;
    }
    check({ displ, strain, stress }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC) {
            return stress !== undefined;
        }
        return displ !== undefined;
    }
    checkCompactionShearBands() {
        if (this.plane[0] && this.plane[1]) {
            // The first and second compaction-shear bands are defined (strike1, dip1 and dipDirection1) and (strike2, dip2 and dipDirection2)
            // if (Number.isNaN(this.params[0].rake) && Number.isNaN(this.params[0].striationTrend)) {
            //     this.striation1 = false
            // }
            // if (Number.isNaN(this.params2.rake) && Number.isNaN(this.params2.striationTrend)) {
            //     this.striation2 = false
            // }
            // if (!this.striation1 && !this.striation2) {
            if (true) {
                // Case 1: general case
                // The striations are not defined for compaction-shear bands 1 and 2
                // In principle we have already checked that the two compaction-shear bands are different
                // Calculate the unit vector parellel to Sigma 2, which is perpendicular to nPlane1 and nPlane2:
                this.nSigma2_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nPlane[0], V: this.nPlane[0] });
                // Calculate the two normalized stress axes (Sigma 1 and Sigma 3) that bisect the angles between the compaction-shear bands
                // angle_nPlane1_nPlane2 in interval (0,PI)
                const angle_nPlane1_nPlane2 = Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: this.nPlane[0], V: this.nPlane[1] }));
                if (Math.abs(angle_nPlane1_nPlane2 - Math.PI / 2) > this.EPS) {
                    // The compaction-shear bands are not perpendicular to each other
                    if (angle_nPlane1_nPlane2 > Math.PI / 2) {
                        // The angle between the 2 normals > PI/2:
                        // Sigma 1 and Sigma 3 are located in the plane generated by nPlane1 and nPlane2 (normal unit vectors point upward)
                        // In principle Sigma 3 bisects the acute angle (< 90°) between the compaction-shear bands
                        // The bisecting line nSigma3_Sc is defined by the sum of normal vectors nPlane1 + nPlane2
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: this.nPlane[0], V: this.nPlane[1] });
                        // note that nSigma3_Sc is always located in the compressional quadrant of the outward hemisphere relative to each of the fault planes
                        // i.e., the scalar product nPlane1 . nSigma3_Sc > 0
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(this.nSigma3_Sc);
                        // nSigma1_Sc = nSigma3_Sc x nSigma2_Sc
                        // The right-handed reference system is defined according to systems S' and S'' (sigma 1, sigma 3, sigma 2)
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nSigma3_Sc, V: this.nSigma2_Sc });
                    }
                    else {
                        // The angle between the 2 normals < PI/2:
                        // In principle Sigma 1 bisects the obtuse angle (> 90°) between the compaction-shear bands 
                        // The bisecting line nSigma1_Sc is defined by the sum of normal vectors nPlane1 + nPlane2
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: this.nPlane[0], V: this.nPlane[1] });
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(this.nSigma1_Sc);
                        // nSigma3_Sc = nSigma2_Sc x nSigma1_Sc
                        // The right-handed reference system is defined according to systems S' and S'' (sigma 1, sigma 3, sigma 2)
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nSigma2_Sc, V: this.nSigma1_Sc });
                    }
                    // ****** let (coordinates1.phi, coordinates1.theta) and (coordinates2.phi, coordinates2.theta) be the spherical coords
                    // of conjugate plaes 1 and 2 in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
                    // This requires using method faultSphericalCoords in class fault
                    // Check that the sense of mouvement is consistent with the orientation of stress axes***
                    if (this.params[0].sensOfMovement !== 'UKN') {
                        // The sense of movement is defined for compaction-shear band 1 (UKN = unknown)
                        // Check consitency of movement 
                        this.csb[0].conjugatePlaneCheckMouvement({
                            noPlane: this.params[0].noPlane,
                            nPlane: this.nPlane[0],
                            coordinates: this.csb[0].fault.sphericalCoords,
                            nStriation: this.csb[0].fault.striation,
                            sensOfMovement: this.params[0].sens_of_movement,
                            nSigma3_Sc: this.nSigma3_Sc,
                            nSigma2_Sc: this.nSigma2_Sc
                        });
                    }
                    if (this.params[1].sensOfMovement !== 'UKN') {
                        // The sense of movement is defined for compaction-shear band 2
                        // Check consitency of movement
                        this.csb[1].conjugatePlaneCheckMouvement({
                            noPlane: this.params[1].noPlane,
                            nPlane: this.nPlane[1],
                            coordinates: this.csb[1].fault.sphericalCoords,
                            nStriation: this.csb[1].fault.striation,
                            sensOfMovement: this.params[1].sens_of_movement,
                            nSigma3_Sc: this.nSigma3_Sc,
                            nSigma2_Sc: this.nSigma2_Sc
                        });
                    }
                }
                else {
                    // The angle between the 2 normals is approximately equal to PI/2:
                    // In this situation, the orientation of Sigma 1 and Sigma 3 can be permuted.
                    // The sense of mouvement of at least one compaction-shear band must be known in order to define the orientation of the stress axes
                    if (this.params[0].sensOfMovement !== 'UKN' || this.params[1].sensOfMovement !== 'UKN') {
                        // Find orientations of Sigma 1, Sigma 2 and Sigma 3, and check for consistency of mouvement if both movements are known.
                    }
                    else {
                        const msg = 'conjugate faults ' + this.params[0].noPlane + ' and ' + this.params[1].noPlane + ' are perpendicular. Please indicate type of movement for at least one plane';
                        throw new Error(msg);
                    }
                }
                // We define 2 orthonormal right-handed reference systems:
                //      S =  (X, Y, Z ) is the geographic reference frame oriented in (East, North, Up) directions.
                //      Sc = (Xc,Yc,Zc) is the principal stress reference frame, parallel to the stress axes obtained from the compaction-shear bands (sigma1_Sc, sigma3_Sc, sigma2_Sc);
                // Calculate transformation matrix Crot from reference system S to Sc, such that:
                //      Vc = Crot V
                // where V and Vc are corresponding vectors in reference systems S and Sc
                // Rotation tensor Crot is defined such that:
                // lines 1, 2, and 3 are defined by unit vectors nSigma1_Sc, nSigma3_Sc, and nSigma2_Sc, in agreement with reference system Sc = (Xc,Yc,Zc)
                this.Crot = [
                    [this.nSigma1_Sc[0], this.nSigma1_Sc[1], this.nSigma1_Sc[2]],
                    [this.nSigma3_Sc[0], this.nSigma3_Sc[1], this.nSigma3_Sc[2]],
                    [this.nSigma2_Sc[0], this.nSigma2_Sc[1], this.nSigma2_Sc[2]]
                ];
            }
            else {}
        }
        else if (this.plane[0] || this.plane[1]) {
            // One of the two compaction-shear bands is defined while the other is not
            // In principle the second plane can be calculated using the striation of the first plane and a line passing by the second plane
            // However, this case is rather unusual.
        }
        else {
            // The two compaction-shear bands are not defined
            throw new Error('conjugate faults ' + this.params[0].noPlane + ' and ' + this.params[1].noPlane + ' are not defined (strike, dip and dip direction');
        }
    }
    cost({ displ, strain, stress, Mrot }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC) {
            // The cost function uses the rotation tensor Crot from reference system S to Sc, calculated in method checkCompactionShearBands
            // The cost function for two compaction-shear bands is defined as the minimum angular rotation between system Sc and the stress tensor in system S' or S'':
            //  S   =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
            //  S'  =  (X', Y', Z' ) is the principal reference frame chosen by the user in the interactive search phase.
            //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
            // Rotation tensors Rrot and RTrot between systems S and S' are such that:
            //  V  = RTrot V'        (RTrot is tensor Rrot transposed)
            //  V' = Rrot  V
            // Rotation tensors Wrot and WTrot between systems S and S'' satisfy : WTrot = RTrot DTrot, such that:
            //  V   = WTrot V''        (WTrot is tensor Wrot transposed)
            //  V'' = Wrot  V
            // The cost method implements a rotation tensor termed Mrot which is equivalent to Rrot or Wrot depending on the calling functions:
            //      Mrot = Rrot in the interactive search phase using integral curves
            //      Mrot = Wrot in the inverse method search phase using Montecarlo (for example)
            // The rotation tensor CFrot between systems Sc and S' (or S'') is such that: Vc = CFrot . V' (or Vc = CFrot . V''), 
            // where CFrot = Crot . MTrot (MTrot = Mrot transposed):
            const CFrot = (0,_types__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: this.Crot, B: (0,_types__WEBPACK_IMPORTED_MODULE_0__.transposeTensor)(Mrot) });
            // const CFrot = multiplyTensors({A: Crot, B: WTrot})
            // The angle of rotation associated to tensor CFrot is defined by the trace tr(CFrot), according to the relation:
            //      tr(CFrot) = 1 + 2 cos(theta)
            // 4 possible right-handed reference systems are considered for CFrot in order to calculate the minimum rotation angle
            return (0,_types_math__WEBPACK_IMPORTED_MODULE_5__.minRotAngleRotationTensor)(CFrot);
        }
    }
}


/***/ }),

/***/ "./lib/data/ConjugateFaults_Kin.ts":
/*!*****************************************!*\
  !*** ./lib/data/ConjugateFaults_Kin.ts ***!
  \*****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "ConjugateFaultsKin": () => (/* binding */ ConjugateFaultsKin)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _Data__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Data */ "./lib/data/Data.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./types */ "./lib/data/types.ts");
/* harmony import */ var _utils_ConjugateFaults__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/ConjugateFaults */ "./lib/utils/ConjugateFaults.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils */ "./lib/utils/index.ts");
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");






/**
 conjugate planes:
 A pair of conjugate planes is defined by two planes whose plane of movement is perpendicular to the intersection line between the planes.
 The plane of movement is defined by the two normal vectors to the fault planes.
 We make the following hypotheses concerning principal stress orientations:
    The compressional axis Sigma 1 is located in the plane of movement and bisects the acute angle (<= 90°) between planes
    The extensional axis Sigma 3 is located in the plane of movement and bisects the obtuse angle (>= 90°) between planes

 Conjugate fault planes are defined in the input file in two consecutive lines.
 Each line specifies all the available data for each conjugate fault plane.

 Several data sets defining two conjugate planes are considered:
 1) Case 1: The geometry and kinematics of the conjugate planes are defined, yet the rake angles are NOT defined
    The orientation of the principal axes are calculated from the geometry of the conjugate fault planes.
    The intermediate axis Sigma 2 is parallel to the intersection line between the conjugate planes;
    The intermediate axis Sigma 2 is perpendicular to the plane of mouvement;
    a) Each plane is defined by a set of 3 parameters, as follows:
        Fault strike: clockwise angle measured from the North direction [0, 360)
        Fault dip: [0, 90]
        Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
    b) The rake angles defining the slip vectors are NOT defined
    c) The sense of mouvement is indicated for each conjugate fault plane :
        For verification purposes, it is recommended to indicate both the dip-slip and strike-slip compoenents, when possible.
          Dip-slip component:
              N = Normal fault,
              I = Inverse fault or thrust
          Strike-slip componenet:
              RL = Right-Lateral fault
              LL = Left-Lateral fault
        Sense of mouvement: N, I, RL, LL, N-RL, N-LL, I-RL, I-LL
2) Case 2: The geometry, the striation (i.e., the rake), and the kinematics of one or both conjugate planes are defined.
        This case is not considered as the striations define redundant information for constraining the stress tensor orientation.
        Conjugate planes with striations can be considered separately as neoformed striated planes for which the friction angle is known

3) Particular case:
    If the striation is known for one of the conjugate faults, then the other plane
    can be defined by one axis that is contained in the plane.
    However, this case would require a special format in the input file, which is inconvenient...

    @category Data
 */
class ConjugateFaultsKin extends _Data__WEBPACK_IMPORTED_MODULE_1__.Data {
    constructor() {
        super(...arguments);
        this.nPlane1 = undefined;
        this.nPlane2 = undefined;
        this.nStriation1 = undefined;
        this.nStriation2 = undefined;
        this.pos = undefined;
        this.problemType = _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC;
        this.strategy = _types__WEBPACK_IMPORTED_MODULE_2__.FractureStrategy.ANGLE;
        this.oriented = true;
        this.EPS = 1e-7;
        this.plane1 = undefined;
        this.plane2 = undefined;
        // Principal directions of the data
        this.nSigma1_Sc = undefined;
        this.nSigma2_Sc = undefined;
        this.nSigma3_Sc = undefined;
        this.Crot = undefined;
        this.cf1 = undefined;
        this.cf2 = undefined;
        this.params1 = undefined;
        this.params2 = undefined;
        this.striation1 = false;
        this.striation2 = false;
    }
    //protected nSigma1_rot: Vector3 = undefined
    //protected nSigma3_rot: Vector3 = undefined
    // params1 and params2 contain data defining conjugate fault planes 1 and 2
    // we have replaced azimuth by strike
    nbLinkedData() {
        return 2;
    }
    initialize(params) {
        let nPlane2Neg;
        if (Number.isNaN(params[0].strike)) {
            throw new Error('Missing strike angle for conjugate fault plane ' + params[0].noPlane);
        }
        if (Number.isNaN(params[0].dip)) {
            throw new Error('Missing dip angle for conjugate fault plane ' + params[0].noPlane);
        }
        if (params[0].dip < 90 && params[0].dipDirection === undefined) {
            throw new Error('Missing dip direction for conjugate fault plane ' + params[0].noPlane);
        }
        if (Number.isNaN(params[1].strike)) {
            throw new Error('Missing strike angle for conjugate fault plane ' + params[1].noPlane);
        }
        if (Number.isNaN(params[1].dip)) {
            throw new Error('Missing dip angle for conjugate fault plane ' + params[1].noPlane);
        }
        if (params[1].dip < 90 && params[1].dipDirection === undefined) {
            throw new Error('Missing dip direction for conjugate fault plane ' + params[1].noPlane);
        }
        // if (this.nPlane1 === this.nPlane2 || this.nPlane1 === constant_x_Vector({k: -1, V: this.nPlane2}) ) {
        //     throw new Error('The two conjugate planes ' + params.noPlane1 + ' and ' + params.noPlane2 + ' are identical')
        // }
        // Check that nPlane and nStriation are unit vectors ***
        this.params1 = {
            noPlane: params[0].noPlane,
            strike: params[0].strike,
            dipDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[0].dipDirection),
            dip: params[0].dip,
            sensOfMovement: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getSensOfMovementFromString)(params[0].typeOfMovement),
            rake: params[0].rake,
            strikeDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[0].strikeDirection)
        };
        this.cf1 = _utils_ConjugateFaults__WEBPACK_IMPORTED_MODULE_3__.ConjugateFaults.create(this.params1);
        // conjugate plane 1 is defined: (strike1, dip1, dipDirection1)
        this.plane1 = true;
        // Calculate the unit vector normal to plane 1: nPlane1
        this.nPlane1 = this.cf1.nPlane;
        // -----------------------------------------
        this.params2 = {
            noPlane: params[1].noPlane,
            strike: params[1].strike,
            dipDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[1].dipDirection),
            dip: params[1].dip,
            sensOfMovement: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getSensOfMovementFromString)(params[1].typeOfMovement),
            rake: params[1].rake,
            strikeDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_4__.getDirectionFromString)(params[1].strikeDirection)
        };
        this.cf2 = _utils_ConjugateFaults__WEBPACK_IMPORTED_MODULE_3__.ConjugateFaults.create(this.params2);
        // conjugate plane 1 is defined: (strike1, dip1, dipDirection1)
        this.plane2 = true;
        // Calculate the unit vector normal to plane 1: nPlane1
        this.nPlane2 = this.cf2.nPlane;
        /** this.nStriation = nStriation
        // this.nPerpStriation = nPerpStriation

        // Check orthogonality
        // const sp = scalarProductUnitVectors({U: nPlane, V: nStriation})
        //*if (Math.abs(sp) >this.EPS) {
            throw new Error(`striation is not on the fault plane. Dot product gives ${sp}`)
        } **/
        return true;
    }
    check({ displ, strain, stress }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC) {
            return stress !== undefined;
        }
        return displ !== undefined;
    }
    checkConjugateFaults() {
        // const f = ConjugatedFaults.create({})
        // return f.
        if (this.plane1 && this.plane2) {
            // The first and second conjugate plane are defined (strike1, dip1 and dipDirection1) and (strike2, dip2 and dipDirection2)
            if (Number.isNaN(this.params1.rake) && Number.isNaN(this.params1.striationTrend)) {
                this.striation1 = false;
            }
            if (Number.isNaN(this.params2.rake) && Number.isNaN(this.params2.striationTrend)) {
                this.striation2 = false;
            }
            if (!this.striation1 && !this.striation2) {
                // Case 1: general case
                // The striations are not defined for the conjugate planes 1 and 2
                // In principle we have already checked that the two conjugate planes are different
                // Calculate the unit vector parellel to Sigma 2, which is perpendicular to nPlane1 and nPlane2:
                this.nSigma2_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nPlane1, V: this.nPlane2 });
                // Calculate the two normalized stress axes (Sigma 1 and Sigma 3) that bisect the angles between the conjugate planes
                // angle_nPlane1_nPlane2 in interval (0,PI)
                const angle_nPlane1_nPlane2 = Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: this.nPlane1, V: this.nPlane2 }));
                if (Math.abs(angle_nPlane1_nPlane2 - Math.PI / 2) > this.EPS) {
                    // The conjugate planes are not perpendicular to each other
                    if (angle_nPlane1_nPlane2 < Math.PI / 2) {
                        // The angle between the 2 normals < PI/2:
                        // Sigma 1 and Sigma 3 are located in the plane generated by nPlane1 and nPlane2  (normal unit vectors point upward)
                        // In principle Sigma 3 bisects the obtuse angle (> 90°) between the conjugate planes
                        // The bisecting line nSigma3_Sc is defined by the sum of normal vectors nPlane1 + nPlane2
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: this.nPlane1, V: this.nPlane2 });
                        // note that nSigma3_Sc is always located in the compressional quadrant of the outward hemisphere relative to each of the fault planes
                        // i.e., the scalar product nPlane1 . nSigma3_Sc > 0
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(this.nSigma3_Sc);
                        // nSigma1_Sc = nSigma3_Sc x nSigma2_Sc
                        // The right-handed reference system is defined according to systems S' and S'' (sigma 1, sigma 3, sigma 2)
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nSigma3_Sc, V: this.nSigma2_Sc });
                    }
                    else {
                        // The angle between the 2 normals > PI/2:
                        // In principle Sigma 1 bisects the acute angle (< 90°) between the conjugate planes 
                        // The bisecting line nSigma1_Sc is defined by the sum of normal vectors nPlane1 + nPlane2
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: this.nPlane1, V: this.nPlane2 });
                        this.nSigma1_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(this.nSigma1_Sc);
                        // nSigma3_Sc = nSigma2_Sc x nSigma1_Sc
                        // The right-handed reference system is defined according to systems S' and S'' (sigma 1, sigma 3, sigma 2)
                        this.nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: this.nSigma2_Sc, V: this.nSigma1_Sc });
                    }
                    // ****** let (coordinates1.phi, coordinates1.theta) and (coordinates2.phi, coordinates2.theta) be the spherical coords
                    // of conjugate plaes 1 and 2 in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
                    // This requires using method faultSphericalCoords in class fault
                    // Check that the sense of mouvement is consistent with the orientation of stress axes***
                    if (this.params1.sensOfMovement !== 'UKN') {
                        // The sense of movement is defined for conjugate plane 1 (UKN = unknown)
                        // Check consitency of movement 
                        this.cf1.conjugatePlaneCheckMouvement({
                            noPlane: this.params1.noPlane,
                            nPlane: this.nPlane1,
                            coordinates: this.cf1.fault.sphericalCoords,
                            nStriation: this.cf1.fault.striation,
                            sensOfMovement: this.params1.sens_of_movement,
                            nSigma3_Sc: this.nSigma3_Sc,
                            nSigma2_Sc: this.nSigma2_Sc
                        });
                    }
                    if (this.params2.sensOfMovement !== 'UKN') {
                        // The sense of movement is defined for conjugate plane 2
                        // Check consitency of movement
                        this.cf2.conjugatePlaneCheckMouvement({
                            noPlane: this.params2.noPlane,
                            nPlane: this.nPlane2,
                            coordinates: this.cf2.fault.sphericalCoords,
                            nStriation: this.cf2.fault.striation,
                            sensOfMovement: this.params2.sens_of_movement,
                            nSigma3_Sc: this.nSigma3_Sc,
                            nSigma2_Sc: this.nSigma2_Sc
                        });
                    }
                }
                else {
                    // The angle between the 2 normals is approximately equal to PI/2:
                    // In this situation, the orientation of Sigma 1 and Sigma 3 can be permuted.
                    // The sense of mouvement of at least one conjugate plane must be known in order to define the orientation of the stress axes
                    if (this.params1.sensOfMovement !== 'UKN' || this.params2.sensOfMovement !== 'UKN') {
                        // Find orientations of Sigma 1, Sigma 2 and Sigma 3, and check for consistency of mouvement if both movements are known.
                    }
                    else {
                        const msg = 'conjugate faults ' + this.params1.noPlane + ' and ' + this.params2.noPlane + ' are perpendicular. Please indicate type of movement for at least one plane';
                        throw new Error(msg);
                    }
                }
                // We define 2 orthonormal right-handed reference systems:
                //      S =  (X, Y, Z ) is the geographic reference frame oriented in (East, North, Up) directions.
                //      Sc = (Xc,Yc,Zc) is the principal stress reference frame, parallel to the stress axes obtained form the conjugate faults (sigma1_Sc, sigma3_Sc, sigma2_Sc);
                // Calculate transformation matrix Crot from reference system S to Sc, such that:
                //      Vc = Crot V
                // where V and Vc are corresponding vectors in reference systems S and Sc
                // Rotation tensor Crot is defined such that:
                // lines 1, 2, and 3 are defined by unit vectors nSigma1_Sc, nSigma3_Sc, and nSigma2_Sc, in agreement with reference system Sc = (Xc,Yc,Zc)
                this.Crot = [
                    [this.nSigma1_Sc[0], this.nSigma1_Sc[1], this.nSigma1_Sc[2]],
                    [this.nSigma3_Sc[0], this.nSigma3_Sc[1], this.nSigma3_Sc[2]],
                    [this.nSigma2_Sc[0], this.nSigma2_Sc[1], this.nSigma2_Sc[2]]
                ];
            }
            else {
                // The striations on one are both conjugate planes are defined (rake, strike direction, and type of movement)
                // 
            }
        }
        else if (this.plane1 || this.plane2) {
            // One of the two conjugate planes is defined while the other is not
            // In principle the second plane can be calculated using the striation of the first plane and a line passing by the second plane
            // However, this case is rather unusual.
        }
        else {
            // The two conjugate planes are not defined
            throw new Error('conjugate faults ' + this.params1.noPlane + ' and ' + this.params2.noPlane + ' are not defined (strike, dip and dip direction');
        }
    }
    cost({ displ, strain, stress, Mrot }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneProblemType.DYNAMIC) {
            // The cost function uses the rotation tensor Crot from reference system S to Sc, calculated in method checkCompactionShearBands
            // The cost function for two conjugate faults is defined as the minimum angular rotation between system Sc and the stress tensor in system S' or S'':
            //  S   =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
            //  S'  =  (X', Y', Z' ) is the principal reference frame chosen by the user in the interactive search phase.
            //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
            // Rotation tensors Rrot and RTrot between systems S and S' are such that:
            //  V  = RTrot V'        (RTrot is tensor Rrot transposed)
            //  V' = Rrot  V
            // Rotation tensors Wrot and WTrot between systems S and S'' satisfy : WTrot = RTrot DTrot, such that:
            //  V   = WTrot V''        (WTrot is tensor Wrot transposed)
            //  V'' = Wrot  V
            // The cost method implements a rotation tensor termed Mrot which is equivalent to Rrot or Wrot depending on the calling functions:
            //      Mrot = Rrot in the interactive search phase using integral curves
            //      Mrot = Wrot in the inverse method search phase using Montecarlo (for example)
            // The rotation tensor CFrot between systems Sc and S' (or S'') is such that: Vc = CFrot . V' (or Vc = CFrot . V''), 
            // where CFrot = Crot . MTrot (MTrot = Mrot transposed):
            const CFrot = (0,_types__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: this.Crot, B: (0,_types__WEBPACK_IMPORTED_MODULE_0__.transposeTensor)(Mrot) });
            // The angle of rotation associated to tensor CFrot is defined by the trace tr(CFrot), according to the relation:
            //      tr(CFrot) = 1 + 2 cos(theta)
            // 4 possible right-handed reference systems are considered for CFrot in order to calculate the minimum rotation angle
            return (0,_types_math__WEBPACK_IMPORTED_MODULE_5__.minRotAngleRotationTensor)(CFrot);
        }
    }
}


/***/ }),

/***/ "./lib/data/Data.ts":
/*!**************************!*\
  !*** ./lib/data/Data.ts ***!
  \**************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "Data": () => (/* binding */ Data)
/* harmony export */ });
/**
 * @brief A Data represents one and only one measure
 * @category Data
 */
class Data {
    constructor() {
        this.weight_ = 1;
        this.active_ = true;
    }
    //private userSpace_ : UserSpace = UserSpace.INVERSE
    weight() {
        return this.weight_;
    }
    set active(a) {
        this.active = a;
    }
    get active() {
        return this.active_;
    }
    description() {
        return undefined;
    }
    setOptions(options) {
        return false;
    }
    nbLinkedData() {
        return 1;
    }
    /**
     * After stress inersion, get the inered data orientation/magnitude/etc for this specific Data
     */
    predict({ displ, strain, stress }) {
        return undefined;
    }
}


/***/ }),

/***/ "./lib/data/DataDescription.ts":
/*!*************************************!*\
  !*** ./lib/data/DataDescription.ts ***!
  \*************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "DataDescription": () => (/* binding */ DataDescription)
/* harmony export */ });
/* harmony import */ var _utils_Fault__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../utils/Fault */ "./lib/utils/Fault.ts");
/* harmony import */ var _Factory__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Factory */ "./lib/data/Factory.ts");


const DataDescription = {
    names: [
        'dataType',
        'strike', 'dip', 'dipDirection',
        'rake', 'strikeDirection', 'striationTrend', 'typeOfMovement',
        'lineTrend', 'linePlunge',
        'deformationPhase', 'relatedWeight',
        'minFrictionAngle', 'maxFrictionAngle', 'minAngleS1n', 'maxAngleS1n' // 12 13 14 15
    ],
    type: [
        (v) => v,
        (v) => parseFloat(v),
        (v) => parseFloat(v),
        (v) => v,
        (v) => parseFloat(v),
        (v) => v,
        (v) => parseFloat(v),
        (v) => v,
        (v) => parseFloat(v),
        (v) => parseFloat(v),
        (v) => parseInt(v),
        (v) => parseFloat(v),
        (v) => parseFloat(v),
        (v) => parseFloat(v),
        (v) => parseFloat(v),
        (v) => parseFloat(v),
    ],
    ranges: [
        (dataType) => _Factory__WEBPACK_IMPORTED_MODULE_1__.DataFactory.exists(dataType),
        (v) => {
            // Strike
            const vv = DataDescription.type[1](v);
            return vv >= 0 && vv < 360;
        },
        (v) => {
            // Dip
            const vv = DataDescription.type[2](v);
            return vv >= 0 && vv <= 90;
        },
        (v) => (0,_utils_Fault__WEBPACK_IMPORTED_MODULE_0__.directionExists)(v),
        (v) => {
            const vv = DataDescription.type[4](v);
            return vv >= 0 && vv <= 90;
        },
        (v) => (0,_utils_Fault__WEBPACK_IMPORTED_MODULE_0__.getDirectionFromString)(v),
        (v) => {
            const vv = DataDescription.type[6](v);
            return vv >= 0 && vv < 360;
        },
        (v) => (0,_utils_Fault__WEBPACK_IMPORTED_MODULE_0__.getSensOfMovementFromString)(v),
        (v) => {
            const vv = DataDescription.type[8](v);
            return vv >= 0 && vv < 360;
        },
        (v) => {
            const vv = DataDescription.type[9](v);
            return vv >= 0 && vv <= 90;
        },
        (v) => {
            const vv = DataDescription.type[10](v);
            return vv > 0;
        },
        (v) => {
            const vv = DataDescription.type[11](v);
            return vv >= 0;
        },
        (v) => {
            const vv = DataDescription.type[12](v);
            return vv >= 0 && vv < 90;
        },
        (v) => {
            const vv = DataDescription.type[13](v);
            return vv >= DataDescription.ranges[12] && vv < 90;
        },
        (v) => {
            const vv = DataDescription.type[14](v);
            return vv > 0 && vv < 90;
        },
        (v) => {
            const vv = DataDescription.type[15](v);
            return vv >= DataDescription.ranges[14] && vv < 90;
        }
    ],
    check({ data, lines, lineNumber }) {
        const result = { status: true, messages: [] };
        const params = [];
        for (let i = 0; i < lines.length; ++i) {
            const line = lines[i];
            const toks = line.split(';').map(s => s.replace(',', '.'));
            if (toks.length !== 16) {
                // TODO: Ben va te faire voir
                throw new Error(`Bad number of columns. Should be 16 and got ${toks.length}`);
            }
            if (this.ranges[0](toks[0]) === false) {
                throw new Error(`Data type named "${toks[0]}" is unknown at line ${lineNumber + i}`);
            }
            const param = {
                dataType: toks[0]
            };
            const desc = data.description();
            desc.mandatory.forEach(index => {
                if (this.ranges[index](toks[index]) === false) {
                    result.status = false;
                    result.messages.push(`for data ${_Factory__WEBPACK_IMPORTED_MODULE_1__.DataFactory.name(data)}, parameter ${this.names[index]} is out of range`);
                }
                else {
                    param[this.names[index]] = this.type[index](toks[index]);
                }
            });
            params.push(param);
        }
        if (result.status === true) {
            // Initialize the data since everything went well
            console.log(params);
            data.initialize(params);
        }
        return result;
    }
};


/***/ }),

/***/ "./lib/data/DataParameters.ts":
/*!************************************!*\
  !*** ./lib/data/DataParameters.ts ***!
  \************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);



/***/ }),

/***/ "./lib/data/ExtensionFracture.ts":
/*!***************************************!*\
  !*** ./lib/data/ExtensionFracture.ts ***!
  \***************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "ExtensionFracture": () => (/* binding */ ExtensionFracture)
/* harmony export */ });
/* harmony import */ var _youwol_math__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @youwol/math */ "@youwol/math");
/* harmony import */ var _youwol_math__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_youwol_math__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _utils_fromAnglesToNormal__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../utils/fromAnglesToNormal */ "./lib/utils/fromAnglesToNormal.ts");
/* harmony import */ var _Data__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Data */ "./lib/data/Data.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./types */ "./lib/data/types.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ../utils */ "./lib/utils/index.ts");






/**
 * @brief Represent an observed and measured joint
 * @category Data
 */
class ExtensionFracture extends _Data__WEBPACK_IMPORTED_MODULE_3__.Data {
    constructor() {
        super(...arguments);
        this.normal = [0, 0, 0];
        this.strategy = _types__WEBPACK_IMPORTED_MODULE_4__.FractureStrategy.ANGLE;
        this.position = undefined;
    }
    description() {
        return {
            mandatory: [1, 2, 3],
            optional: [10, 11]
        };
    }
    initialize(params) {
        if (Number.isNaN(params[0].strike)) {
            throw new Error('Missing azimuth angle for ExtensionFracture');
        }
        if (Number.isNaN(params[0].dip)) {
            throw new Error('Missing dip angle for ExtensionFracture');
        }
        if (params[0].dip < 90 && params[0].dipDirection === undefined) {
            throw new Error('Missing dip direction for ExtensionFracture');
        }
        // Convert into normal
        this.normal = (0,_utils_fromAnglesToNormal__WEBPACK_IMPORTED_MODULE_2__.fromAnglesToNormal)({
            strike: params[0].strike,
            dip: params[0].dip,
            dipDirection: (0,_utils__WEBPACK_IMPORTED_MODULE_5__.getDirectionFromString)(params[0].dipDirection)
        });
        return true;
    }
    check({ displ, strain, stress }) {
        return stress !== undefined;
    }
    // This version does not consider the case in which the stress shape ratio R is close to zero (i.e., Sigma 2 = Sigma 3) 
    //      and any plane containing Sigma 1 is consistent with the stress tensor solution
    cost({ displ, strain, stress }) {
        // [xx, xy, xz, yy, yz, zz]
        const sigma = [stress[0][0], stress[0][1], stress[0][2], stress[1][1], stress[1][2], stress[2][2]];
        // eigen = function calculating the 3 normalized eigenvectors (Sigma_1, Sigma_2, Sigma_3) of the stress tensor ??
        // vectors is formated like: [S1x, S1y, S1z, S2x..., S3z]
        const { values, vectors } = (0,_youwol_math__WEBPACK_IMPORTED_MODULE_0__.eigen)(sigma);
        const s = [vectors[6], vectors[7], vectors[8]];
        // dot = scalar product between 2 unitary vectors Sigma_3 . normal = cos(angle)
        const dot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.scalarProductUnitVectors)({ U: s, V: this.normal });
        switch (this.strategy) {
            case _types__WEBPACK_IMPORTED_MODULE_4__.FractureStrategy.DOT: return 1 - Math.abs(dot);
            // Sigma 3 can be oriented in two opposite directions, thus to calculate the minimum angle we take the dot product as positive.
            default: return Math.acos(Math.abs(dot)) / Math.PI;
        }
    }
}


/***/ }),

/***/ "./lib/data/Factory.ts":
/*!*****************************!*\
  !*** ./lib/data/Factory.ts ***!
  \*****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "DataFactory": () => (/* binding */ DataFactory)
/* harmony export */ });
/* harmony import */ var _ExtensionFracture__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./ExtensionFracture */ "./lib/data/ExtensionFracture.ts");
/* harmony import */ var _StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./StriatedPlane_Kin */ "./lib/data/StriatedPlane_Kin.ts");
/* harmony import */ var _StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./StriatedPlane_Friction1 */ "./lib/data/StriatedPlane_Friction1.ts");
/* harmony import */ var _StriatedPlane_Friction2__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./StriatedPlane_Friction2 */ "./lib/data/StriatedPlane_Friction2.ts");
/* harmony import */ var _CompactionShearBands_Kin__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./CompactionShearBands_Kin */ "./lib/data/CompactionShearBands_Kin.ts");
/* harmony import */ var _ConjugateFaults_Kin__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./ConjugateFaults_Kin */ "./lib/data/ConjugateFaults_Kin.ts");






/* eslint @typescript-eslint/no-explicit-any: off -- need to have any here for the factory */
var DataFactory;
(function (DataFactory) {
    const map_ = new Map();
    DataFactory.bind = (obj, name = '') => {
        name.length === 0 ? map_.set(obj.name, obj) : map_.set(name, obj);
    };
    DataFactory.create = (name, params = undefined) => {
        const M = map_.get(name);
        if (M) {
            return new M(params);
        }
        return undefined;
    };
    DataFactory.exists = (name) => {
        return map_.get(name) !== undefined;
    };
    DataFactory.names = () => {
        return Array.from(map_.keys());
    };
    DataFactory.name = (data) => {
        return data.constructor.name;
    };
})(DataFactory || (DataFactory = {}));
DataFactory.bind(_ExtensionFracture__WEBPACK_IMPORTED_MODULE_0__.ExtensionFracture, 'Extension Fracture');
DataFactory.bind(_StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneKin, 'Striated Plane');
DataFactory.bind(_StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_2__.StriatedPlaneFriction1, 'Striated Plane Friction1');
DataFactory.bind(_StriatedPlane_Friction2__WEBPACK_IMPORTED_MODULE_3__.StriatedPlaneFriction2, 'Striated Plane Friction2');
DataFactory.bind(_CompactionShearBands_Kin__WEBPACK_IMPORTED_MODULE_4__.CompactionShearBandsKin, 'Compaction Shear Bands');
DataFactory.bind(_ConjugateFaults_Kin__WEBPACK_IMPORTED_MODULE_5__.ConjugateFaultsKin, 'Conjugate Faults');


/***/ }),

/***/ "./lib/data/StriatedPlane_Friction1.ts":
/*!*********************************************!*\
  !*** ./lib/data/StriatedPlane_Friction1.ts ***!
  \*********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "StriatedPlaneFriction1": () => (/* binding */ StriatedPlaneFriction1)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./StriatedPlane_Kin */ "./lib/data/StriatedPlane_Kin.ts");


/**
 * @category Data
 */
class StriatedPlaneFriction1 extends _StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneKin {
    constructor() {
        super(...arguments);
        this.cohesionRock_ = 0;
        this.frictionAngleRock_ = 0;
        this.weightFriction_ = 1;
    }
    set cohesionRock(c) {
        this.cohesionRock_ = c;
    }
    set frictionAngleRock(f) {
        this.frictionAngleRock_ = f;
        if (this.frictionAngleRock_ <= this.EPS) {
            // A positive friction angle has to be defined prior to stress tensor inversion
            throw ('For friction analysis choose frictionAngleRock > 0 ');
        }
    }
    // Not used yet!!!
    set weightFriction(w) {
        this.weightFriction_ = w;
    }
    initialize(params) {
        if (super.initialize(params) === false) {
            return false;
        }
        if (params[0].cohesion) {
            this.cohesionRock = params[0].cohesion;
        }
        if (params[0].friction) {
            this.frictionAngleRock = params[0].friction;
        }
        if (params[0].weightFriction) {
            this.weightFriction = params[0].weightFriction;
        }
        return true;
    }
    cost({ displ, strain, stress }) {
        // For each striated fault the misfit distance is defined in terms of an angular distance between the stress vector F
        // and the closest axis Fkf that satisfies both the kinematical and frictional costraints.
        // In other words, Fkf is such that its projection on the fault plane is parallel to the measured striation 
        // and it is located along or above the friction line in the Mohr Circle diagram.
        // The normal stress is calculated by shifting the origin of the normalized Mohr circle toward the left, such that the friction law passes by the new origin
        // This condition allows to calculate friction angles for the total stress vectors that can be directly compared with the rock friction angle
        // Moreover, this condition is consistent with a residual friction law for shear faulting
        // The misfit distance corresponds to an angular distance between the stress vector F and Fkf
        let misfitDistance = Number.MAX_VALUE;
        let stress2;
        let stress_Shifted_Sigma_n;
        let F;
        let Flocal;
        const frictionStriaUnitVecLocal = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        let magStressShiftSigma_n;
        let deltaNormalStress;
        // deltaNormalStress = Shift of the normalized Mohr circle along the normal stress axis
        //      such that the friction line intersects the origin of the plane (normal stress, shear stress)
        //      deltaNormalStress > 0, according to rock mechanics sign convention : Compressional stresses > 0
        deltaNormalStress = this.cohesionRock_ / Math.tan(this.frictionAngleRock_);
        // Let (xl,yl,zl) be a local right-handed reference frame fixed to the fault plane, where:
        //      xl = unit vector pointing toward the striation
        //      yl = unit vector perpendicular to the striation and located in the fault plane
        //      zl = fault normal (i.e., pointing upward)
        // Note that (xl,yl,zl) is defined by 3 orthonormal vectors: (fault.striation, fault.e_perp_striation, fault.normal)
        // frictionStriaUnitVecLocal = constant unit vector (Fkfo) in local reference frame (xl,yl,zl) such that:
        //      1) Its projection on the fault plane is parallel to the measured striation
        //      2) It is located along the friction line
        frictionStriaUnitVecLocal[0] = Math.sin(this.frictionAngleRock_);
        frictionStriaUnitVecLocal[1] = 0;
        // The normal stress is defined according to the rock mechanics sign convention : Compressional stresses > 0
        frictionStriaUnitVecLocal[2] = Math.cos(this.frictionAngleRock_);
        // ------------- FOR (before)
        //==============  Stress analysis using continuum mechanics sign convention : Compressional stresses < 0
        // In principle, principal stresses are negative: (sigma 1, sigma 2, sigma 3) = (-1, -R, 0) 
        // Calculate total stress vector F:  
        stress2 = (0,_types__WEBPACK_IMPORTED_MODULE_0__.tensor_x_Vector)({ T: stress, V: this.nPlane });
        // stress_Shifted_Sigma_n = Stress vector obtained by adding the shift of the normal stress component deltaNormalStress
        //      i.e., in a 'frictional' reference frame such that the Mohr Coulomb line intersects the origin
        stress_Shifted_Sigma_n = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: stress2, V: (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -deltaNormalStress, V: this.nPlane }) });
        magStressShiftSigma_n = (0,_types__WEBPACK_IMPORTED_MODULE_0__.vectorMagnitude)(stress_Shifted_Sigma_n);
        if (magStressShiftSigma_n > this.EPS) { // stressMag > Epsilon ***
            // In principle the stress magnitude > 0 ***
            // F = shifted stress vector (compression < 0), normalized for angular calculations
            F = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(stress_Shifted_Sigma_n, magStressShiftSigma_n);
            //==============  Friction analysis using rock mechanics sign convention : Compressional stresses > 0
            // Calculate the normalized stress vector F in local reference frame (xl,yl,zl)
            Flocal[0] = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nStriation });
            // For angular calculations the stress componenet parallel to yl is taken as positive:
            //  the sign of yl has no influence on the misfit distance
            Flocal[1] = Math.abs((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nPerpStriation }));
            // For angular calculations the normal stress is taken as positive consistently with frictionStriaUnitVecLocal[2] > 0
            // This tranformation is equivalent to a reflexion of the stress vector relative to the plane (xl,yl)
            // In other words, the sign of the normal stress is inverted from negative to positive
            Flocal[2] = Math.abs((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nPlane }));
            // Calculate the colatitude angle theta_z, in interval [0, PI/2] since yl>= 0 
            //      yl = cos(theta_z)
            let theta_z = Math.acos(Flocal[1]);
            if (Math.abs(theta_z) > this.EPS) {
                // Calculate the azimuthal angle phi_z necessary for determining the misfit distance in terms of the orientation of the stress vector
                // phi_z is in interval [-PI/2, PI/2] depending on the sign of xl
                //      xl = sin(theta_z) * sin(phi_z)
                let phi_z = Math.asin(Flocal[0] / Math.sin(theta_z));
                if (phi_z >= this.frictionAngleRock_) {
                    // Case 1: The stress vector satisfies the frictional criteriun
                    // The angular distance to the closest axis Fkf is located along a great circle passing by yl.
                    // Note that the misfit angular difference is lower than the angular difference between measured and calculated striae
                    misfitDistance = (Math.PI / 2) - theta_z;
                }
                else if (phi_z > 0) {
                    // Case 2: 
                    // xl > 0, thus the angular difference between measured and calculated striation < PI/2
                    // Note that the stress vector may satisfy (or not) the frictional criteriun;
                    // The misfit is defined by the angular deviation of the stress relative to frictionStriaUnitVecLocal.
                    // More precisely, frictionStriaUnitVecLocal is the closest axis satisfying the kinematic and frictional criteria
                    misfitDistance = Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: Flocal, V: frictionStriaUnitVecLocal }));
                }
                else {
                    // Case 3: xl <= 0, thus the angular difference between measured and calculated striation >= PI/2
                    // In principle, these faults should be eliminated from the solution set by assigning a misfit greater than PI/2
                    misfitDistance = Math.max(Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: Flocal, V: frictionStriaUnitVecLocal })), Math.PI / 2);
                }
            }
            else {
                // stressUnitVecLocal is parallel to yl: stressUnitVecLocal[1] = 1
                misfitDistance = Math.PI / 2;
            }
        }
        else {
            // The friction angle is zero and the plane is perpendicular to Sigma_3
            // In such situation, the fault should be eliminated from the solution set by assigning a misfit of PI/2
            misfitDistance = Math.PI / 2;
        }
        return misfitDistance;
    }
}


/***/ }),

/***/ "./lib/data/StriatedPlane_Friction2.ts":
/*!*********************************************!*\
  !*** ./lib/data/StriatedPlane_Friction2.ts ***!
  \*********************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "StriatedPlaneFriction2": () => (/* binding */ StriatedPlaneFriction2)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./StriatedPlane_Friction1 */ "./lib/data/StriatedPlane_Friction1.ts");


/**
 * @category Data
 */
class StriatedPlaneFriction2 extends _StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneFriction1 {
    cost({ displ, strain, stress }) {
        // For each striated fault the misfit distance is defined in terms of an angular distance between the stress vector F
        // and the closest axis Fkf that satisfies both the kinematical and frictional costraints.
        // In other words, Fkf is such that its projection on the fault plane is parallel to the measured striation 
        // and it is located along or above the friction line in the Mohr Circle diagram.
        // The normal stress is calculated by shifting the origin of the normalized Mohr circle toward the left, such that the friction law passes by the new origin
        // This condition allows to calculate friction angles for the total stress vectors that can be directly compared with the rock friction angle
        // Moreover, this condition is consistent with a residual friction law for shear faulting
        // The misfit distance corresponds to an angular distance between the stress vector F and Fkf
        let misfitDistance = Number.MAX_VALUE;
        let stress2;
        let stress_Shifted_Sigma_n;
        let F;
        let Flocal;
        const frictionStriaUnitVecLocal = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        let magStressShiftSigma_n;
        let deltaNormalStress;
        if (this.frictionAngleRock_ > this.EPS) {
            // deltaNormalStress = Shift of the normalized Mohr circle along the normal stress axis
            //      such that the friction line intersects the origin of the plane (normal stress, shear stress)
            //      deltaNormalStress > 0, according to rock mechanics sign convention : Compressional stresses > 0
            deltaNormalStress = this.cohesionRock_ / Math.tan(this.frictionAngleRock_);
        }
        else {
            // A positive friction angle has to be defined prior to stress tensor inversion
            throw ('For friction analysis choose frictionAngleRock > 0 ');
        }
        // Let (xl,yl,zl) be a local right-handed reference frame fixed to the fault plane, where:
        //      xl = unit vector pointing toward the striation
        //      yl = unit vector perpendicular to the striation and located in the fault plane
        //      zl = fault normal (i.e., pointing upward)
        // Note that (xl,yl,zl) is defined by 3 orthonormal vectors: (fault.striation, fault.e_perp_striation, fault.normal)
        // frictionStriaUnitVecLocal = constant unit vector (Fkfo) in local reference frame (xl,yl,zl) such that:
        //      1) Its projection on the fault plane is parallel to the measured striation
        //      2) It is located along the friction line
        frictionStriaUnitVecLocal[0] = Math.sin(this.frictionAngleRock_);
        frictionStriaUnitVecLocal[1] = 0;
        // The normal stress is defined according to the rock mechanics sign convention : Compressional stresses > 0
        frictionStriaUnitVecLocal[2] = Math.cos(this.frictionAngleRock_);
        //==============  Stress analysis using continuum mechanics sign convention : Compressional stresses < 0
        // In principle, principal stresses are negative: (sigma 1, sigma 2, sigma 3) = (-1, -R, 0) 
        // Calculate total stress vector F:  
        stress2 = (0,_types__WEBPACK_IMPORTED_MODULE_0__.tensor_x_Vector)({ T: stress, V: this.nPlane });
        // stress_Shifted_Sigma_n = Stress vector obtained by adding the shift of the normal stress component deltaNormalStress
        //      i.e., in a 'frictional' reference frame such that the Mohr Coulomb line intersects the origin
        stress_Shifted_Sigma_n = (0,_types__WEBPACK_IMPORTED_MODULE_0__.add_Vectors)({ U: stress2, V: (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -deltaNormalStress, V: this.nPlane }) });
        magStressShiftSigma_n = (0,_types__WEBPACK_IMPORTED_MODULE_0__.vectorMagnitude)(stress_Shifted_Sigma_n);
        if (magStressShiftSigma_n > this.EPS) { // stressMag > Epsilon ***
            // In principle the stress magnitude > 0 ***
            // F = shifted stress vector (compression < 0), normalized for angular calculations
            F = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(stress_Shifted_Sigma_n, magStressShiftSigma_n);
            //==============  Friction analysis using rock mechanics sign convention : Compressional stresses > 0
            // Calculate the normalized stress vector F in local reference frame (xl,yl,zl)
            Flocal[0] = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nStriation });
            // For angular calculations the stress componenet parallel to yl is taken as positive:
            //  the sign of yl has no influence on the misfit distance
            Flocal[1] = Math.abs((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nPerpStriation }));
            // For angular calculations the normal stress is taken as positive consistently with frictionStriaUnitVecLocal[2] > 0
            // This tranformation is equivalent to a reflexion of the stress vector relative to the plane (xl,yl)
            // In other words, the sign of the normal stress is inverted from negative to positive
            Flocal[2] = Math.abs((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: F, V: this.nPlane }));
            // Calculate the colatitude angle theta_z, in interval [0, PI/2] since yl>= 0 
            //      yl = cos(theta_z)
            let theta_z = Math.acos(Flocal[1]);
            if (Math.abs(theta_z) > this.EPS) {
                // Calculate the azimuthal angle phi_z necessary for determining the misfit distance in terms of the orientation of the stress vector
                // phi_z is in interval [-PI/2, PI/2] depending on the sign of xl
                //      xl = sin(theta_z) * sin(phi_z)
                let phi_z = Math.asin(Flocal[0] / Math.sin(theta_z));
                if (phi_z >= this.frictionAngleRock_) {
                    // Case 1: The stress vector satisfies the frictional criteriun
                    // The angular distance to the closest axis Fkf is located along a great circle passing by yl.
                    // Note that the misfit angular difference is lower than the angular difference between measured and calculated striae
                    misfitDistance = (Math.PI / 2) - theta_z;
                }
                else if (phi_z > 0) {
                    // Case 2: 
                    // xl > 0, thus the angular difference between measured and calculated striation < PI/2
                    // Note that the stress vector may satisfy (or not) the frictional criteriun;
                    // The misfit is defined by the angular deviation of the stress relative to frictionStriaUnitVecLocal.
                    // More precisely, frictionStriaUnitVecLocal is the closest axis satisfying the kinematic and frictional criteria
                    misfitDistance = Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: Flocal, V: frictionStriaUnitVecLocal }));
                }
                else {
                    // Case 3: xl <= 0, thus the angular difference between measured and calculated striation >= PI/2
                    // In principle, these faults should be eliminated from the solution set by assigning a misfit greater than PI/2
                    misfitDistance = Math.max(Math.acos((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: Flocal, V: frictionStriaUnitVecLocal })), Math.PI / 2);
                }
            }
            else {
                // stressUnitVecLocal is parallel to yl: stressUnitVecLocal[1] = 1
                misfitDistance = Math.PI / 2;
            }
        }
        else {
            // The friction angle is zero and the plane is perpendicular to Sigma_3
            // In such situation, the fault should be eliminated from the solution set by assigning a misfit of PI/2
            misfitDistance = Math.PI / 2;
        }
        return misfitDistance;
    }
}


/***/ }),

/***/ "./lib/data/StriatedPlane_Kin.ts":
/*!***************************************!*\
  !*** ./lib/data/StriatedPlane_Kin.ts ***!
  \***************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "StriatedPlaneKin": () => (/* binding */ StriatedPlaneKin)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _Data__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Data */ "./lib/data/Data.ts");
/* harmony import */ var _types_mechanics__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../types/mechanics */ "./lib/types/mechanics.ts");
/* harmony import */ var _utils_Fault__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/Fault */ "./lib/utils/Fault.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./types */ "./lib/data/types.ts");





/**
 * @category Data
 */
class StriatedPlaneKin extends _Data__WEBPACK_IMPORTED_MODULE_1__.Data {
    constructor() {
        super(...arguments);
        this.nPlane = undefined;
        this.nStriation = undefined;
        this.pos = undefined;
        this.problemType = _types__WEBPACK_IMPORTED_MODULE_4__.StriatedPlaneProblemType.DYNAMIC;
        this.strategy = _types__WEBPACK_IMPORTED_MODULE_4__.FractureStrategy.ANGLE;
        this.oriented = true;
        this.EPS = 1e-7;
    }
    initialize(params) {
        if (Number.isNaN(params[0].strike)) {
            throw new Error('Missing strike angle for StriatedPlaneKin');
        }
        if (Number.isNaN(params[0].dip)) {
            throw new Error('Missing dip angle for StriatedPlaneKin');
        }
        if (params[0].dip < 90 && Number.isNaN(params[0].dipDirection)) {
            throw new Error('Missing dip direction for StriatedPlaneKin');
        }
        // Check that nPlane and nStriation are unit vectors
        const { nPlane, nStriation, nPerpStriation } = _utils_Fault__WEBPACK_IMPORTED_MODULE_3__.Fault.create({
            strike: params[0].strike,
            dipDirection: this.getMapDirection(params[0].dipDirection),
            dip: params[0].dip,
            sensOfMovement: this.getSensOfMovement(params[0].typeOfMovement),
            rake: params[0].rake,
            strikeDirection: this.getMapDirection(params[0].strikeDirection)
        });
        this.nPlane = nPlane;
        this.nStriation = nStriation;
        this.nPerpStriation = nPerpStriation;
        // Check orthogonality
        const sp = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nPlane, V: nStriation });
        if (Math.abs(sp) > this.EPS) {
            throw new Error(`striation is not on the fault plane. Dot product gives ${sp}`);
        }
        return true;
    }
    check({ displ, strain, stress }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_4__.StriatedPlaneProblemType.DYNAMIC) {
            return stress !== undefined;
        }
        return displ !== undefined;
    }
    cost({ displ, strain, stress }) {
        if (this.problemType === _types__WEBPACK_IMPORTED_MODULE_4__.StriatedPlaneProblemType.DYNAMIC) {
            // For the first implementation, use the W&B hyp.
            // let d = tensor_x_Vector({T: stress, V: this.nPlane}) // Cauchy
            // d = normalizeVector(d)
            // Calculate shear stress parameters
            // Calculate the magnitude of the shear stress vector in reference system S
            const { shearStress, normalStress, shearStressMag } = (0,_types_mechanics__WEBPACK_IMPORTED_MODULE_2__.faultStressComponents)({ stressTensor: stress, normal: this.nPlane });
            let cosAngularDifStriae = 0;
            if (shearStressMag > 0) { // shearStressMag > Epsilon would be more realistic ***
                // nShearStress = unit vector parallel to the shear stress (i.e. representing the calculated striation)
                let nShearStress = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizeVector)(shearStress, shearStressMag);
                // The angular difference is calculated using the scalar product: 
                // nShearStress . nStriation = |nShearStress| |nStriation| cos(angularDifStriae) = 1 . 1 . cos(angularDifStriae)
                // cosAngularDifStriae = cos(angular difference between calculated and measured striae)
                cosAngularDifStriae = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nShearStress, V: this.nStriation });
            }
            else {
                // The calculated shear stress is zero (i.e., the fault plane is parallel to a principal stress)
                // In such situation we may consider that the calculated striation can have any direction.
                // Nevertheless, the plane should not display striations as the shear stress is zero.
                // Thus, in principle the plane is not compatible with the stress tensor, and it should be eliminated from the analysis
                // In suchh case, the angular difference is taken as PI
                cosAngularDifStriae = -1;
            }
            if (this.strategy === _types__WEBPACK_IMPORTED_MODULE_4__.FractureStrategy.ANGLE) {
                // The misfit is defined by the angular difference (in radians) between measured and calculated striae
                if (this.oriented) {
                    // The sense of the striation is known
                    return Math.acos(cosAngularDifStriae);
                }
                else {
                    // The sense of the striation is not known. Thus, we choose the sens that minimizes the angular difference 
                    // and is more compatible with the observed striation.
                    return Math.acos(Math.abs(cosAngularDifStriae));
                }
            }
            else {
                // The misfit is defined by the the cosine of the angular difference between measured and calculated striae
                if (this.oriented) {
                    return 0.5 - cosAngularDifStriae / 2;
                }
                else {
                    return 0.5 - Math.abs(cosAngularDifStriae) / 2;
                }
            }
        }
        throw new Error('Kinematic not yet available');
    }
    getMapDirection(s) {
        if (!mapDirection.has(s)) {
            throw new Error(`Direction ${s} is not defined (or incorrectly defined)`);
        }
        return mapDirection.get(s);
    }
    getSensOfMovement(s) {
        if (!mapSensOfMovement.has(s)) {
            throw new Error(`Sens of movement ${s} is not defined (or incorrectly defined)`);
        }
        return mapSensOfMovement.get(s);
    }
}
// ----------------------------------------------------
const mapDirection = new Map();
mapDirection.set("E", 0 /* Direction.E */);
mapDirection.set("N", 2 /* Direction.N */);
mapDirection.set("NE", 4 /* Direction.NE */);
mapDirection.set("NW", 7 /* Direction.NW */);
mapDirection.set("S", 3 /* Direction.S */);
mapDirection.set("SE", 5 /* Direction.SE */);
mapDirection.set("SW", 6 /* Direction.SW */);
mapDirection.set("W", 1 /* Direction.W */);
const mapSensOfMovement = new Map();
mapSensOfMovement.set("Inverse", 2 /* SensOfMovement.I */);
mapSensOfMovement.set("Inverse - Left Lateral", 8 /* SensOfMovement.I_LL */);
mapSensOfMovement.set("Inverse - Right Lateral", 7 /* SensOfMovement.I_RL */);
mapSensOfMovement.set("Left Lateral", 4 /* SensOfMovement.LL */);
mapSensOfMovement.set("Normal", 1 /* SensOfMovement.N */);
mapSensOfMovement.set("Normal - Left Lateral", 6 /* SensOfMovement.N_LL */);
mapSensOfMovement.set("Normal - Right Lateral", 5 /* SensOfMovement.N_RL */);
mapSensOfMovement.set("Right Lateral", 3 /* SensOfMovement.RL */);
mapSensOfMovement.set("Unknown", 9 /* SensOfMovement.UKN */);


/***/ }),

/***/ "./lib/data/index.ts":
/*!***************************!*\
  !*** ./lib/data/index.ts ***!
  \***************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "CompactionShearBandsKin": () => (/* reexport safe */ _CompactionShearBands_Kin__WEBPACK_IMPORTED_MODULE_9__.CompactionShearBandsKin),
/* harmony export */   "ConjugateFaultsKin": () => (/* reexport safe */ _ConjugateFaults_Kin__WEBPACK_IMPORTED_MODULE_10__.ConjugateFaultsKin),
/* harmony export */   "Data": () => (/* reexport safe */ _Data__WEBPACK_IMPORTED_MODULE_0__.Data),
/* harmony export */   "DataDescription": () => (/* reexport safe */ _DataDescription__WEBPACK_IMPORTED_MODULE_1__.DataDescription),
/* harmony export */   "DataFactory": () => (/* reexport safe */ _Factory__WEBPACK_IMPORTED_MODULE_3__.DataFactory),
/* harmony export */   "ExtensionFracture": () => (/* reexport safe */ _ExtensionFracture__WEBPACK_IMPORTED_MODULE_5__.ExtensionFracture),
/* harmony export */   "FractureStrategy": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_4__.FractureStrategy),
/* harmony export */   "StriatedPlaneFriction1": () => (/* reexport safe */ _StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_7__.StriatedPlaneFriction1),
/* harmony export */   "StriatedPlaneFriction2": () => (/* reexport safe */ _StriatedPlane_Friction2__WEBPACK_IMPORTED_MODULE_8__.StriatedPlaneFriction2),
/* harmony export */   "StriatedPlaneKin": () => (/* reexport safe */ _StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_6__.StriatedPlaneKin),
/* harmony export */   "StriatedPlaneProblemType": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_4__.StriatedPlaneProblemType)
/* harmony export */ });
/* harmony import */ var _Data__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Data */ "./lib/data/Data.ts");
/* harmony import */ var _DataDescription__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./DataDescription */ "./lib/data/DataDescription.ts");
/* harmony import */ var _DataParameters__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./DataParameters */ "./lib/data/DataParameters.ts");
/* harmony import */ var _Factory__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Factory */ "./lib/data/Factory.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./types */ "./lib/data/types.ts");
/* harmony import */ var _ExtensionFracture__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./ExtensionFracture */ "./lib/data/ExtensionFracture.ts");
/* harmony import */ var _StriatedPlane_Kin__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./StriatedPlane_Kin */ "./lib/data/StriatedPlane_Kin.ts");
/* harmony import */ var _StriatedPlane_Friction1__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./StriatedPlane_Friction1 */ "./lib/data/StriatedPlane_Friction1.ts");
/* harmony import */ var _StriatedPlane_Friction2__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./StriatedPlane_Friction2 */ "./lib/data/StriatedPlane_Friction2.ts");
/* harmony import */ var _CompactionShearBands_Kin__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./CompactionShearBands_Kin */ "./lib/data/CompactionShearBands_Kin.ts");
/* harmony import */ var _ConjugateFaults_Kin__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./ConjugateFaults_Kin */ "./lib/data/ConjugateFaults_Kin.ts");











// export * from './StyloliteInterface'
// export * from './StyloliteTeeth'
// export * from './SecondaryFault'
// export * from './FocalMechanism_Kin'


/***/ }),

/***/ "./lib/data/types.ts":
/*!***************************!*\
  !*** ./lib/data/types.ts ***!
  \***************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "FractureStrategy": () => (/* binding */ FractureStrategy),
/* harmony export */   "StriatedPlaneProblemType": () => (/* binding */ StriatedPlaneProblemType)
/* harmony export */ });
/**
 * @category Data
 */
var FractureStrategy;
(function (FractureStrategy) {
    FractureStrategy[FractureStrategy["ANGLE"] = 0] = "ANGLE";
    FractureStrategy[FractureStrategy["DOT"] = 1] = "DOT";
})(FractureStrategy || (FractureStrategy = {}));
/**
 * - DYNAMIC is related to forces (or stresses)
 * - KINEMATIC is related to displacement field
 * @category Data
 */
var StriatedPlaneProblemType;
(function (StriatedPlaneProblemType) {
    StriatedPlaneProblemType[StriatedPlaneProblemType["DYNAMIC"] = 0] = "DYNAMIC";
    StriatedPlaneProblemType[StriatedPlaneProblemType["KINEMATIC"] = 1] = "KINEMATIC";
})(StriatedPlaneProblemType || (StriatedPlaneProblemType = {}));


/***/ }),

/***/ "./lib/index.ts":
/*!**********************!*\
  !*** ./lib/index.ts ***!
  \**********************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "CompactionShearBandsKin": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.CompactionShearBandsKin),
/* harmony export */   "ConjugateFaultsKin": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.ConjugateFaultsKin),
/* harmony export */   "Curve3D": () => (/* reexport safe */ _analysis__WEBPACK_IMPORTED_MODULE_0__.Curve3D),
/* harmony export */   "Data": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.Data),
/* harmony export */   "DataDescription": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.DataDescription),
/* harmony export */   "DataFactory": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.DataFactory),
/* harmony export */   "Direction": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.Direction),
/* harmony export */   "EquipotentialCurve": () => (/* reexport safe */ _analysis__WEBPACK_IMPORTED_MODULE_0__.EquipotentialCurve),
/* harmony export */   "ExtensionFracture": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.ExtensionFracture),
/* harmony export */   "Fault": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.Fault),
/* harmony export */   "FractureStrategy": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.FractureStrategy),
/* harmony export */   "GridSearch": () => (/* reexport safe */ _search__WEBPACK_IMPORTED_MODULE_4__.GridSearch),
/* harmony export */   "IntegralCurve": () => (/* reexport safe */ _analysis__WEBPACK_IMPORTED_MODULE_0__.IntegralCurve),
/* harmony export */   "InverseMethod": () => (/* reexport safe */ _InverseMethod__WEBPACK_IMPORTED_MODULE_6__.InverseMethod),
/* harmony export */   "MasterStress": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.MasterStress),
/* harmony export */   "MohrCoulombCurve": () => (/* reexport safe */ _analysis__WEBPACK_IMPORTED_MODULE_0__.MohrCoulombCurve),
/* harmony export */   "MonteCarlo": () => (/* reexport safe */ _search__WEBPACK_IMPORTED_MODULE_4__.MonteCarlo),
/* harmony export */   "PoleCoords": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.PoleCoords),
/* harmony export */   "SearchMethodFactory": () => (/* reexport safe */ _search__WEBPACK_IMPORTED_MODULE_4__.SearchMethodFactory),
/* harmony export */   "SensOfMovement": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.SensOfMovement),
/* harmony export */   "SphericalCoords": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.SphericalCoords),
/* harmony export */   "StressTensor": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.StressTensor),
/* harmony export */   "StriatedPlaneFriction1": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneFriction1),
/* harmony export */   "StriatedPlaneFriction2": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneFriction2),
/* harmony export */   "StriatedPlaneKin": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneKin),
/* harmony export */   "StriatedPlaneProblemType": () => (/* reexport safe */ _data__WEBPACK_IMPORTED_MODULE_1__.StriatedPlaneProblemType),
/* harmony export */   "add_Vectors": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.add_Vectors),
/* harmony export */   "angularDifStriations": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.angularDifStriations),
/* harmony export */   "arcCircle": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.arcCircle),
/* harmony export */   "cloneMatrix3x3": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.cloneMatrix3x3),
/* harmony export */   "cloneMisfitCriteriunSolution": () => (/* reexport safe */ _InverseMethod__WEBPACK_IMPORTED_MODULE_6__.cloneMisfitCriteriunSolution),
/* harmony export */   "constant_x_Vector": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.constant_x_Vector),
/* harmony export */   "createDefaultSolution": () => (/* reexport safe */ _InverseMethod__WEBPACK_IMPORTED_MODULE_6__.createDefaultSolution),
/* harmony export */   "crossProduct": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.crossProduct),
/* harmony export */   "decodeCSV": () => (/* reexport safe */ _io__WEBPACK_IMPORTED_MODULE_5__.decodeCSV),
/* harmony export */   "deg2rad": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.deg2rad),
/* harmony export */   "directionExists": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.directionExists),
/* harmony export */   "faultParams": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.faultParams),
/* harmony export */   "faultStressComponents": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.faultStressComponents),
/* harmony export */   "getDirectionFromString": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.getDirectionFromString),
/* harmony export */   "getSensOfMovementFromString": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.getSensOfMovementFromString),
/* harmony export */   "lineSphericalCoords": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.lineSphericalCoords),
/* harmony export */   "minRotAngleRotationTensor": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.minRotAngleRotationTensor),
/* harmony export */   "mohrCircleLine": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.mohrCircleLine),
/* harmony export */   "mohrCirclePoint": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.mohrCirclePoint),
/* harmony export */   "multiplyTensors": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.multiplyTensors),
/* harmony export */   "newMatrix3x3": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.newMatrix3x3),
/* harmony export */   "newMatrix3x3Identity": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.newMatrix3x3Identity),
/* harmony export */   "newPoint3D": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.newPoint3D),
/* harmony export */   "newVector3D": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.newVector3D),
/* harmony export */   "normalVector": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.normalVector),
/* harmony export */   "normalizeVector": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.normalizeVector),
/* harmony export */   "normalizedCrossProduct": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.normalizedCrossProduct),
/* harmony export */   "properRotationTensor": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.properRotationTensor),
/* harmony export */   "rad2deg": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.rad2deg),
/* harmony export */   "rotationParamsFromRotTensor": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.rotationParamsFromRotTensor),
/* harmony export */   "scalarProduct": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.scalarProduct),
/* harmony export */   "scalarProductUnitVectors": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.scalarProductUnitVectors),
/* harmony export */   "setValueInUnitInterval": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.setValueInUnitInterval),
/* harmony export */   "signedAngularDifStriations": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.signedAngularDifStriations),
/* harmony export */   "spherical2unitVectorCartesian": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.spherical2unitVectorCartesian),
/* harmony export */   "stressTensorDelta": () => (/* reexport safe */ _search__WEBPACK_IMPORTED_MODULE_4__.stressTensorDelta),
/* harmony export */   "stressTensorPrincipalAxes": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.stressTensorPrincipalAxes),
/* harmony export */   "tensor_x_Vector": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.tensor_x_Vector),
/* harmony export */   "transposeTensor": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.transposeTensor),
/* harmony export */   "trend2phi": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.trend2phi),
/* harmony export */   "trimAll": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_3__.trimAll),
/* harmony export */   "unitVectorCartesian2Spherical": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.unitVectorCartesian2Spherical),
/* harmony export */   "vectorMagnitude": () => (/* reexport safe */ _types__WEBPACK_IMPORTED_MODULE_2__.vectorMagnitude)
/* harmony export */ });
/* harmony import */ var _analysis__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./analysis */ "./lib/analysis/index.ts");
/* harmony import */ var _data__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./data */ "./lib/data/index.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./types */ "./lib/types/index.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./utils */ "./lib/utils/index.ts");
/* harmony import */ var _search__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./search */ "./lib/search/index.ts");
/* harmony import */ var _io__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./io */ "./lib/io/index.ts");
/* harmony import */ var _InverseMethod__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./InverseMethod */ "./lib/InverseMethod.ts");









/***/ }),

/***/ "./lib/io/decodeCSV.ts":
/*!*****************************!*\
  !*** ./lib/io/decodeCSV.ts ***!
  \*****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "decodeCSV": () => (/* binding */ decodeCSV)
/* harmony export */ });
/* harmony import */ var _data__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../data */ "./lib/data/index.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../utils */ "./lib/utils/index.ts");


function decodeCSV(lines) {
    const datas = [];
    let dataType = '';
    for (let i = 0; i < lines.length; ++i) {
        if (i === 1) {
            continue;
        }
        const line = (0,_utils__WEBPACK_IMPORTED_MODULE_1__.trimAll)(lines[i].trim());
        if (line.length === 0) {
            continue;
        }
        const r = line.split(';').filter(v => v.length !== 0).map(s => s.replace(',', '.'));
        if (r.length === 0) {
            continue;
        }
        if (i === 0) {
            dataType = r[0];
            console.log(`Data type is "${dataType}"`);
            continue;
        }
        if (r.length === 6) {
            const n = [parseFloat(r[0]), parseFloat(r[1]), parseFloat(r[2])];
            const s = [parseFloat(r[3]), parseFloat(r[4]), parseFloat(r[5])];
            const data = _data__WEBPACK_IMPORTED_MODULE_0__.DataFactory.create(dataType, { nPlane: n, nStriation: s });
            datas.push(data);
        }
    }
    return datas;
}


/***/ }),

/***/ "./lib/io/index.ts":
/*!*************************!*\
  !*** ./lib/io/index.ts ***!
  \*************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "decodeCSV": () => (/* reexport safe */ _decodeCSV__WEBPACK_IMPORTED_MODULE_0__.decodeCSV)
/* harmony export */ });
/* harmony import */ var _decodeCSV__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./decodeCSV */ "./lib/io/decodeCSV.ts");

// export * from './decodeCSV_Angles'


/***/ }),

/***/ "./lib/search/DebugSearch.ts":
/*!***********************************!*\
  !*** ./lib/search/DebugSearch.ts ***!
  \***********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "DebugSearch": () => (/* binding */ DebugSearch)
/* harmony export */ });
/* harmony import */ var _InverseMethod__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../InverseMethod */ "./lib/InverseMethod.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");


class DebugSearch {
    setInteractiveSolution({ rot, stressRatio }) {
    }
    run(data, misfitCriteriaSolution) {
        const newSolution = (0,_InverseMethod__WEBPACK_IMPORTED_MODULE_0__.cloneMisfitCriteriunSolution)(misfitCriteriaSolution);
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                const stress = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
                const misfit = data.reduce((previous, current) => previous + current.cost({ stress }), 0) / data.length;
                if (misfit < newSolution.misfit) {
                    newSolution.misfit = misfit;
                    newSolution.rotationMatrixD = stress;
                    newSolution.stressRatio = undefined;
                    newSolution.stressTensorSolution = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)();
                }
            }
        }
        return newSolution;
    }
}


/***/ }),

/***/ "./lib/search/Factory.ts":
/*!*******************************!*\
  !*** ./lib/search/Factory.ts ***!
  \*******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "SearchMethodFactory": () => (/* binding */ SearchMethodFactory)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _DebugSearch__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./DebugSearch */ "./lib/search/DebugSearch.ts");
/* harmony import */ var _GridSearch__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./GridSearch */ "./lib/search/GridSearch.ts");
/* harmony import */ var _MonteCarlo__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./MonteCarlo */ "./lib/search/MonteCarlo.ts");




var SearchMethodFactory;
(function (SearchMethodFactory) {
    const map_ = new Map();
    SearchMethodFactory.bind = (obj, name = '') => {
        name.length === 0 ? map_.set(obj.name, obj) : map_.set(name, obj);
    };
    SearchMethodFactory.create = (name, params = undefined) => {
        const M = map_.get(name);
        if (M) {
            const searchMethod = new M(params);
            // to be filled
            const ist = params.interactiveStressTensor;
            const st = new _types__WEBPACK_IMPORTED_MODULE_0__.StressTensor({
                trendS1: ist.trendS1,
                trendS3: ist.trendS3,
                plungeS1: ist.plungeS1,
                plungeS3: ist.plungeS3,
                masterStress: ist.masterStress === 'Sigma1' ? _types__WEBPACK_IMPORTED_MODULE_0__.MasterStress.Sigma1 : _types__WEBPACK_IMPORTED_MODULE_0__.MasterStress.Sigma3,
                stressRatio: ist.stressRatio
            });
            searchMethod.setInteractiveSolution({ rot: st.Rrot, stressRatio: st.stressRatio });
            return searchMethod;
        }
        return undefined;
    };
    SearchMethodFactory.exists = (name) => {
        return map_.get(name) !== undefined;
    };
    SearchMethodFactory.names = () => {
        return Array.from(map_.keys());
    };
})(SearchMethodFactory || (SearchMethodFactory = {}));
SearchMethodFactory.bind(_GridSearch__WEBPACK_IMPORTED_MODULE_2__.GridSearch, 'Grid Search');
SearchMethodFactory.bind(_DebugSearch__WEBPACK_IMPORTED_MODULE_1__.DebugSearch, 'Debug Search');
SearchMethodFactory.bind(_MonteCarlo__WEBPACK_IMPORTED_MODULE_3__.MonteCarlo, 'Monte Carlo');


/***/ }),

/***/ "./lib/search/GridSearch.ts":
/*!**********************************!*\
  !*** ./lib/search/GridSearch.ts ***!
  \**********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "GridSearch": () => (/* binding */ GridSearch)
/* harmony export */ });
/* harmony import */ var _InverseMethod__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../InverseMethod */ "./lib/InverseMethod.ts");
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "./lib/search/utils.ts");




/**
 * @category Search-Method
 */
class GridSearch {
    constructor({ deltaGridAngle = 1, GridAngleHalfIntervalS = 30, stressRatioHalfInterval = 0.2, deltaStressRatio = 0.01,
    // Rrot=newMatrix3x3Identity(),
    // stressRatio=0.5
     } = {}) {
        this.deltaGridAngle = 0;
        this.GridAngleHalfIntervalS = 0;
        this.stressRatioHalfInterval = 0;
        this.deltaStressRatio = 0;
        this.Rrot = undefined;
        this.stressRatio0 = 0;
        this.deltaGridAngle = deltaGridAngle;
        this.GridAngleHalfIntervalS = GridAngleHalfIntervalS;
        this.stressRatioHalfInterval = stressRatioHalfInterval;
        this.deltaStressRatio = deltaStressRatio;
        // this.Rrot = Rrot
        // this.stressRatio0 = stressRatio
    }
    setInteractiveSolution({ rot, stressRatio }) {
        this.Rrot = rot;
        this.stressRatio0 = stressRatio;
    }
    /**
     * Example:
     * ```ts
     * const searchMethod = new GridSearch()
     * searchMethod.setInteractiveSolution({Rrot: rot, stressRatio: 0.5})
     *
     * const criterion = new Etchecopar({faultSet, maxNbFault: 100})
     *
     * const initSolution = createDefaultSolution(criterion)
     *
     * for (let i=0; i<100; ++i) {
     *      const solution = searchMethod.run(initSolution)
     * }
     * ```
     */
    run(data, misfitCriteriaSolution) {
        // The optimum stress tensor is calculated by exploring the stress orientations and the stress ratio around the approximate solution S0
        // obtained by the user during the interactive analysis of flow lines on the sphere, Mohr circle diagram, and histogram of signed angular deviations.
        // More precisely, the minimization function is calculated in the nodes of a four-dimmensional grid that sweeps the area around S0
        // The angular node interval englobes the angular interval around the estimated stress directions defined by the user
        let nodesAngleInterval = Math.ceil(this.GridAngleHalfIntervalS / this.deltaGridAngle);
        // The stress ratio node interval englobes the stress ratio interval around the estimated value defined by the user
        let nodesStressRatioInterval = Math.ceil(this.stressRatioHalfInterval / this.deltaStressRatio);
        let deltaGridAngleRad = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.deg2rad)(this.deltaGridAngle);
        let DTrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let Drot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let WTrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let Wrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let inc = 0;
        const newSolution = (0,_InverseMethod__WEBPACK_IMPORTED_MODULE_0__.cloneMisfitCriteriunSolution)(misfitCriteriaSolution);
        for (let i = -nodesAngleInterval; i <= nodesAngleInterval; i++) {
            // Angular variation around axis X': ROLL
            let deltaPhi = i * deltaGridAngleRad;
            let cosDeltaPhi = Math.cos(deltaPhi);
            let sinDeltaPhi = Math.sin(deltaPhi);
            for (let j = -nodesAngleInterval; j <= nodesAngleInterval; j++) {
                // Angular variation around axis Y': PITCH
                let deltaTheta = j * deltaGridAngleRad;
                let cosDeltaTheta = Math.cos(deltaTheta);
                let sinDeltaTheta = Math.sin(deltaTheta);
                for (let k = -nodesAngleInterval; k <= nodesAngleInterval; k++) {
                    // Angular variation around axis Z': YAW
                    const deltaAlpha = k * deltaGridAngleRad;
                    let cosDeltaAlpha = Math.cos(deltaAlpha);
                    let sinDeltaAlpha = Math.sin(deltaAlpha);
                    // Calculate rotation tensors Drot and DTrot between systems S' and S'' such that:
                    //  V'  = DTrot V''        (DTrot is tensor Drot transposed)
                    //  V'' = Drot  V'
                    DTrot = rotationTensorDT(cosDeltaPhi, sinDeltaPhi, cosDeltaTheta, sinDeltaTheta, cosDeltaAlpha, sinDeltaAlpha);
                    Drot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(DTrot);
                    // To analyse the distribution of rotation axes in space: (these instructions can be deleted once the analysis is done)
                    // The cartesian and spherical coords of a unit vector corresponding to the rotation axis are determined 
                    // from the components of the tensor definning a proper rotation
                    // let {rotAxis, rotAxisSpheCoords, rotMag} = rotationParamsFromRotTensor({rotTensor: DTrot}) // **
                    // Calculate rotation tensors Wrot and WTrot between systems S and S'': WTrot = RTrot DTrot, such that:
                    //  V   = WTrot V''        (WTrot is tensor Wrot transposed)
                    //  V'' = Wrot  V
                    //  S   =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
                    //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
                    WTrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.multiplyTensors)({ A: (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(this.Rrot), B: DTrot });
                    //  Wrot = Drot Rrot
                    Wrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(WTrot);
                    for (let l = -nodesStressRatioInterval; l <= nodesStressRatioInterval; l++) {
                        // Stress ratio variation around R = (S2-S3)/(S1-S3)
                        let stressRatio = this.stressRatio0 + l * this.deltaStressRatio;
                        if (stressRatio >= 0 && stressRatio <= 1) { // The strees ratio is in interval [0,1]
                            // Calculate the stress tensor STdelta in reference frame S from the stress tensor in reference frame S''
                            let STdelta = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.stressTensorDelta)(stressRatio, Wrot, WTrot);
                            const misfit = data.reduce((previous, current) => {
                                //console.log(current.cost({stress: STdelta}))
                                return previous + current.cost({ stress: STdelta });
                            }, 0) / data.length;
                            if (misfit < newSolution.misfit) {
                                newSolution.misfit = misfit;
                                newSolution.rotationMatrixD = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(Drot);
                                newSolution.rotationMatrixW = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(Wrot);
                                newSolution.stressRatio = stressRatio;
                                newSolution.stressTensorSolution = STdelta;
                            }
                            inc++;
                        }
                    }
                }
            }
        }
        return newSolution;
        // To analyse the rotation axis for the best solution: 
        // The cartesian and spherical coords of a unit vector corresponding to the rotation axis are determined 
        // from the components of the tensor definning a proper rotation
        // let {rotAxis, rotAxisSpheCoords, rotMag} = rotationParamsFromRotTensor(DTrot) // **
    }
}
// --------------- Hidden to users
function rotationTensorDT(cosDeltaPhi, sinDeltaPhi, cosDeltaTheta, sinDeltaTheta, cosDeltaAlpha, sinDeltaAlpha) {
    // Calculate the rotation tensor DT between reference frame S' and S'', such that:
    //  V'  = DT V''        (DT is tensor D transposed)
    //  V'' = D  V'
    //  S' = (X',Y',Z') is the principal stress reference frame obtained by the user from the interactive analysis, parallel to (sigma_1, sigma_3, sigma_2);
    //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
    // The columns of matrix D are given by the unit vectors parallel to X1'', X2'', and X3'' defined in reference system S':
    const DT = (0,_types_math__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
    // Sigma_1 axis: Unit vector e1''
    DT[0][0] = cosDeltaPhi * cosDeltaTheta;
    DT[1][0] = sinDeltaPhi * cosDeltaTheta;
    DT[2][0] = -sinDeltaTheta;
    // Sigma_3 axis: Unit vector e2''
    DT[0][1] = -sinDeltaPhi * cosDeltaAlpha + cosDeltaPhi * sinDeltaTheta * sinDeltaAlpha;
    DT[1][1] = cosDeltaPhi * cosDeltaAlpha + sinDeltaPhi * sinDeltaTheta * sinDeltaAlpha;
    DT[2][1] = cosDeltaTheta * sinDeltaAlpha;
    // Sigma_2 axis: Unit vector e3''
    DT[0][2] = sinDeltaPhi * sinDeltaAlpha + cosDeltaPhi * sinDeltaTheta * cosDeltaAlpha;
    DT[1][2] = -cosDeltaPhi * sinDeltaAlpha + sinDeltaPhi * sinDeltaTheta * cosDeltaAlpha;
    DT[2][2] = cosDeltaTheta * cosDeltaAlpha;
    return DT;
}


/***/ }),

/***/ "./lib/search/MonteCarlo.ts":
/*!**********************************!*\
  !*** ./lib/search/MonteCarlo.ts ***!
  \**********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "MonteCarlo": () => (/* binding */ MonteCarlo)
/* harmony export */ });
/* harmony import */ var _InverseMethod__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../InverseMethod */ "./lib/InverseMethod.ts");
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "./lib/search/utils.ts");



class MonteCarlo {
    setInteractiveSolution({ rot, stressRatio }) {
        this.Rrot = rot;
        this.stressRatio0 = stressRatio;
    }
    constructor({ stressRatio = 0.5, stressRatioHalfInterval = 0.25, rotAngleHalfInterval = 0.1, nbRandomTrials = 100, Rrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3Identity)() } = {}) {
        this.Rrot = undefined;
        this.RTrot = undefined;
        this.rotAngleHalfInterval = rotAngleHalfInterval;
        this.nbRandomTrials = nbRandomTrials;
        this.stressRatio0 = stressRatio;
        this.stressRatioHalfInterval = stressRatioHalfInterval;
        this.Rrot = Rrot;
        this.RTrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(this.Rrot);
    }
    run(data, misfitCriteriaSolution) {
        // The optimum stress tensor is calculated by exploring the stress orientations and the stress ratio around the approximate solution S0
        // obtained by the user during the interactive analysis of flow lines on the sphere, Mohr circle diagram, and histogram of signed angular deviations.
        // More precisely, the minimization function is calculated for a set of stress tensors whose orientations are rotated around axes 
        // defined by a Montecarlo algorithm distributed "quasi-homogeneously" on the sphere surface.
        // The magnitude of rotations and the stress ratio are also defined by random variables calculated within specified intervals.
        // Stress ratio variation stressRatioHalfInterval around R = (S2-S3)/(S1-S3), is constrained to interval [0,1]
        let stressRatioMin = Math.max(0, Math.abs(this.stressRatio0) - this.stressRatioHalfInterval);
        let stressRatioMax = Math.min(1, Math.abs(this.stressRatio0) + this.stressRatioHalfInterval);
        let stressRatioEffectiveInterval = stressRatioMax - stressRatioMin;
        // console.log(this.stressRatio0, stressRatioMin, stressRatioMax, stressRatioEffectiveInterval)
        let DTrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let Drot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let WTrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let Wrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.newMatrix3x3)();
        let rotAxisSpheCoords = new _types__WEBPACK_IMPORTED_MODULE_1__.SphericalCoords;
        let inc = 0;
        console.log('Starting the montecarlo search...');
        let changed = false;
        const newSolution = (0,_InverseMethod__WEBPACK_IMPORTED_MODULE_0__.cloneMisfitCriteriunSolution)(misfitCriteriaSolution);
        for (let i = 0; i <= this.nbRandomTrials; i++) {
            // For each trial, a rotation axis in the unit sphere is calculated from a uniform random distribution.
            // phi = random variable representing azimuth [0, 2PI)
            rotAxisSpheCoords.phi = Math.random() * 2 * Math.PI;
            // theta = random variable representing the colatitude [0, PI)
            //      the arcos function ensures a uniform distribution for theta from a random value:
            rotAxisSpheCoords.theta = Math.acos(2 * Math.random() - 1);
            let rotAxis = (0,_types__WEBPACK_IMPORTED_MODULE_1__.spherical2unitVectorCartesian)(rotAxisSpheCoords);
            // We only consider positive rotation angles around each rotation axis, since the whole sphere is covered by angles (phi,theta)
            let rotAngle = Math.random() * this.rotAngleHalfInterval;
            // Calculate rotation tensors Drot and DTrot between systems S' and S'' such that:
            //  V'  = DTrot V''        (DTrot is tensor Drot transposed)
            //  V'' = Drot  V'
            DTrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.properRotationTensor)({ nRot: rotAxis, angle: rotAngle });
            Drot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(DTrot);
            // Calculate rotation tensors Wrot and WTrot between systems S and S'': WTrot = RTrot DTrot, such that:
            //  V   = WTrot V''        (WTrot is tensor Wrot transposed)
            //  V'' = Wrot  V
            //  S   =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
            //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
            WTrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.multiplyTensors)({ A: this.RTrot, B: DTrot });
            //  Wrot = Drot Rrot
            Wrot = (0,_types__WEBPACK_IMPORTED_MODULE_1__.transposeTensor)(WTrot);
            // Stress ratio variation around R = (S2-S3)/(S1-S3)
            let stressRatio = stressRatioMin + Math.random() * stressRatioEffectiveInterval; // The strees ratio is in interval [0,1]
            // Calculate the stress tensor STdelta in reference frame S from the stress tensor in reference frame S''
            // STdelta is defined according to the continuum mechanics sign convention : compression < 0
            let STdelta = (0,_utils__WEBPACK_IMPORTED_MODULE_2__.stressTensorDelta)(stressRatio, Wrot, WTrot);
            const misfit = data.reduce((previous, current) => {
                return previous + current.cost({ stress: STdelta, rot: Wrot });
            }, 0) / data.length;
            if (misfit < newSolution.misfit) {
                newSolution.misfit = misfit;
                newSolution.rotationMatrixD = (0,_types__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(Drot);
                newSolution.rotationMatrixW = (0,_types__WEBPACK_IMPORTED_MODULE_1__.cloneMatrix3x3)(Wrot);
                newSolution.stressRatio = stressRatio;
                newSolution.stressTensorSolution = STdelta;
            }
            // const misfitSum  = misfitCriteriaSolution.criterion.value(STdelta)
            // if (misfitSum < misfitCriteriaSolution.misfitSum) {
            //     misfitCriteriaSolution.misfitSum      = misfitSum
            //     misfitCriteriaSolution.rotationMatrixD = cloneMatrix3x3(Drot)
            //     misfitCriteriaSolution.rotationMatrixW = cloneMatrix3x3(Wrot)
            //     misfitCriteriaSolution.stressRatio    = stressRatio
            //     changed = true
            // }
            inc++;
        }
        return newSolution;
    }
}
/*
private monteCarloSearch() {
        // The optimum stress tensor is calculated by exploring the stress orientations and the stress ratio around the approximate solution S0
        // obtained by the user during the interactive analysis of flow lines on the sphere, Mohr circle diagram, and histogram of signed angular deviations.
        // More precisely, the minimization function is calculated for a set of stress tensors whose orientations are rotated around axes
        // defined by a Montecarlo algorithm distributed "quasi-homogeneously" on the sphere surface.
        // The magnitude of rotations and the stress ratio are also defined by random variables calculated within specified intervals.

        // Stress ratio variation stressRatioHalfInterval around R = (S2-S3)/(S1-S3), is constrained to interval [0,1]
        let stressRatioMin = Math.max(0, this.stressRatio0 - this.stressRatioHalfInterval )
        let stressRatioMax = Math.min(1, this.stressRatio0 + this.stressRatioHalfInterval )

        let stressRatioEffectiveInterval = stressRatioMax - stressRatioMin
         
        let DTrot: Matrix3x3   = newMatrix3x3()
        let Drot:  Matrix3x3   = newMatrix3x3()
        let WTrot: Matrix3x3   = newMatrix3x3()
        let Wrot:  Matrix3x3   = newMatrix3x3()

        let rotAxisSpheCoords: SphericalCoords

        let inc = 0

        console.log('Starting the grid search...')
        
        for (let i = 0; i <= this.nRandomTrials; i++) {
            // For each trial, a rotation axis in the unit sphere is calculated from a uniform random distribution.

            // phi = random variable representing azimuth [0, 2PI)
            rotAxisSpheCoords.phi = Math.random() * 2 * Math.PI
            // theta = random variable representing the colatitude [0, PI)
            //      the arcos function ensures a uniform distribution for theta from a random value:
            rotAxisSpheCoords.theta = Math.acos( Math.random() * 2 * Math.PI)

            let rotAxis = spherical2unitVectorCartesian(rotAxisSpheCoords)

            // We only consider positive rotation angles around each rotation axis, since the whole sphere is covered by angles (phi,theta)
            let rotAngle = Math.random() * this.rotAngleHalfInterval
                
            // Calculate rotation tensors Drot and DTrot between systems S' and S'' such that:
            //  V'  = DTrot V''        (DTrot is tensor Drot transposed)
            //  V'' = Drot  V'
            DTrot = properRotationTensor(rotAxis, rotAngle)
            Drot  = transposeTensor(DTrot)
            // Calculate rotation tensors Wrot and WTrot between systems S and S'': WTrot = RTrot DT, such that:
            //  V   = WTrot V''        (WTrot is tensor Wrot transposed)
            //  V'' = Wrot  V
            //  S   =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
            //  S'' =  (X'', Y'', Z'' ) is the principal reference frame for a fixed node in the search grid (sigma_1, sigma_3, sigma_2)
            WTrot = multiplyTensors({A: this.RTrot, B: DTrot })
            Wrot  = transposeTensor( WTrot )

            // Stress ratio variation around R = (S2-S3)/(S1-S3)
            let stressRatio = stressRatioMin + Math.random() * stressRatioEffectiveInterval // The strees ratio is in interval [0,1]
            // Calculate the stress tensor STdelta in reference frame S from the stress tensor in reference frame S''
            let STdelta = stressTensorDelta(stressRatio, Wrot, WTrot)

            this.misfitCriteriaSolution.forEach( bestSolution => {
                const misfitSum  = bestSolution.criterion({
                    rotationMatrixDrot: Drot,
                    rotationMatrixWrot: Wrot,
                    stressRatio: stressRatio,
                    faultSet: this.faultSet
                }) // number
                if (misfitSum < bestSolution.misfitSum) {
                    bestSolution.misfitSum      = misfitSum
                    bestSolution.rotationMatrixDrot = cloneMatrix3x3(Drot)
                    bestSolution.rotationMatrixWrot = cloneMatrix3x3(Wrot)
                    bestSolution.stressRatio    = stressRatio
                    console.log('cur solution:', [i,j,k,l])
                }
            })

            // if (inc%100 === 0) {
                // console.log(inc, [i,j,k,l])
                console.log(inc, rotAxisSpheCoords)
            // }

            inc++
        }
        // To analyse the rotation axis for the best solution:
        // The cartesian and spherical coords of a unit vector corresponding to the rotation axis are determined
        // from the components of the tensor definning a proper rotation
        // let {rotAxis, rotAxisSpheCoords, rotMag} = rotationParamsFromRotTensor(DTrot) // **
    }
*/


/***/ }),

/***/ "./lib/search/SearchMethod.ts":
/*!************************************!*\
  !*** ./lib/search/SearchMethod.ts ***!
  \************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);



/***/ }),

/***/ "./lib/search/index.ts":
/*!*****************************!*\
  !*** ./lib/search/index.ts ***!
  \*****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "GridSearch": () => (/* reexport safe */ _GridSearch__WEBPACK_IMPORTED_MODULE_3__.GridSearch),
/* harmony export */   "MonteCarlo": () => (/* reexport safe */ _MonteCarlo__WEBPACK_IMPORTED_MODULE_4__.MonteCarlo),
/* harmony export */   "SearchMethodFactory": () => (/* reexport safe */ _Factory__WEBPACK_IMPORTED_MODULE_1__.SearchMethodFactory),
/* harmony export */   "stressTensorDelta": () => (/* reexport safe */ _utils__WEBPACK_IMPORTED_MODULE_2__.stressTensorDelta)
/* harmony export */ });
/* harmony import */ var _SearchMethod__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./SearchMethod */ "./lib/search/SearchMethod.ts");
/* harmony import */ var _Factory__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Factory */ "./lib/search/Factory.ts");
/* harmony import */ var _utils__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./utils */ "./lib/search/utils.ts");
/* harmony import */ var _GridSearch__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./GridSearch */ "./lib/search/GridSearch.ts");
/* harmony import */ var _MonteCarlo__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./MonteCarlo */ "./lib/search/MonteCarlo.ts");







/***/ }),

/***/ "./lib/search/utils.ts":
/*!*****************************!*\
  !*** ./lib/search/utils.ts ***!
  \*****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "stressTensorDelta": () => (/* binding */ stressTensorDelta)
/* harmony export */ });
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");

/**
 * @category Search-Method
 * @param stressRatio
 * @param Wrot
 * @param WTrot
 * @returns
 */
function stressTensorDelta(stressRatio, Wrot, WTrot) {
    // Calculate the stress tensor STdelta in reference frame S from the stress tensor in reference frame S'':
    //      STdelta = WTrot STPdelta Wrot
    //
    // where
    //
    //      S'' = (X'',Y'',Z'') is the principal stress reference frame in the search grid node(i,j,k,l), parallel to (sigma_1, sigma_3, sigma_2);
    //      S   = (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
    //      STPdelta = Stress Tensor in the Principal stress reference frame (i.e. diagonal tensor with eigenvalues (1,0,StressRatio).
    //      The principal stress values are negative since stress calculations are done using the continuum mechanics convention.
    const sigma = [-1, 0, -stressRatio];
    let STPdelta = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.stressTensorPrincipalAxes)(sigma);
    // T1 = WTrot STPdelta; this tensor multiplication can be optimized since STPdelta is diagonal with eigenvalues (-1, 0, -StressRatio).
    let T1 = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
    T1 = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: WTrot, B: STPdelta });
    // STdelta = T1 Wrot = WTrot STPdelta Wrot
    let STdelta = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: T1, B: Wrot });
    return STdelta;
    // for (let m = 0; m < numberStressInversions; m++) {
    //     // The user may stipulate 1 or 2 different stress inversion methods for the same fault set
    //     // This option allows to compare inversion solutions using different misfit criteria
    //     switch(MisfitCriteriun[m]){
    //         case 1: {
    //             // Angular deviation (Etchecopar et al. 1981)
    //             // FirstNode allows to initialize the mimimisation function with its value for the fist node in the grid
    //             misfitAngularDeviation( firstNode, m )
    //         }
    //         case 2: { 
    //             //Minimum angle of rotation of the tensor that brings the shear stress parallel to the striation (Gephart & Forsyth 1984)
    //             misfitMinimumAngleTensorRotation( firstNode, m )
    //         }
    //     }                                
    // }
}


/***/ }),

/***/ "./lib/types/GenericCurve.ts":
/*!***********************************!*\
  !*** ./lib/types/GenericCurve.ts ***!
  \***********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);



/***/ }),

/***/ "./lib/types/MohrPoint.ts":
/*!********************************!*\
  !*** ./lib/types/MohrPoint.ts ***!
  \********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);



/***/ }),

/***/ "./lib/types/PoleCoords.ts":
/*!*********************************!*\
  !*** ./lib/types/PoleCoords.ts ***!
  \*********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "PoleCoords": () => (/* binding */ PoleCoords)
/* harmony export */ });
class PoleCoords {
    constructor({ trend = 0, plunge = 0 } = {}) {
        this.trend_ = 0;
        this.plunge_ = 0;
        // check bounds of theta and phi if any
        this.trend_ = trend ? trend : 0;
        this.plunge_ = plunge ? plunge : 0;
    }
    get trend() {
        return this.trend_;
    }
    set trend(v) {
        // check bounds of theta
        this.trend_ = v;
    }
    get plunge() {
        return this.plunge_;
    }
    set plunge(v) {
        this.plunge_ = v;
    }
}


/***/ }),

/***/ "./lib/types/SphericalCoords.ts":
/*!**************************************!*\
  !*** ./lib/types/SphericalCoords.ts ***!
  \**************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "SphericalCoords": () => (/* binding */ SphericalCoords)
/* harmony export */ });
class SphericalCoords {
    constructor({ phi = 0, theta = 0 } = {}) {
        this.theta_ = 0;
        this.phi_ = 0;
        // check bounds of theta and phi if any
        this.theta_ = theta ? theta : 0;
        this.phi_ = phi ? phi : 0;
    }
    /**
     * @see constructor
     */
    static create({ phi = 0, theta = 0 } = {}) {
        return new SphericalCoords({ phi, theta });
    }
    get theta() {
        return this.theta_;
    }
    set theta(v) {
        // check bounds of theta
        this.theta_ = v;
    }
    get phi() {
        return this.phi_;
    }
    set phi(v) {
        this.phi_ = v;
    }
}


/***/ }),

/***/ "./lib/types/StressTensor.ts":
/*!***********************************!*\
  !*** ./lib/types/StressTensor.ts ***!
  \***********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "MasterStress": () => (/* binding */ MasterStress),
/* harmony export */   "StressTensor": () => (/* binding */ StressTensor)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");

var MasterStress;
(function (MasterStress) {
    MasterStress[MasterStress["Sigma1"] = 0] = "Sigma1";
    MasterStress[MasterStress["Sigma3"] = 1] = "Sigma3";
})(MasterStress || (MasterStress = {}));
class StressTensor {
    constructor({ trendS1, plungeS1, trendS3, plungeS3, masterStress, stressRatio }) {
        // ========================================================
        this.phiS1 = 0;
        this.phiS3 = 0;
        this.thetaS1 = 0;
        this.thetaS3 = 0;
        /**
         * The new stress tensor is defined in a new reference sytem S' = (X',Y',Z') = (sigma 1, sigma 3, sigma 2)
         * Sigma 1 and sigma 3 axes can be defined interactively in the sphere (prefered solution) or chosen in a predefined range.
         *
         * The stress axis Sigma_1 is defined by two angles in PoleCoords: trend and plunge.
         *      trend  = azimuth of sigma 1 in interval [0, 360), measured clockwise from the North
         *      plunge = vertical angle between the horizontal plane and the sigma 1 axis (positive downward), in interval [-90,90]
         * @param param0
         */
        this.poleS1 = new _types__WEBPACK_IMPORTED_MODULE_0__.PoleCoords();
        /**
         * The stress axis Sigma_3 is defined by two angles in PoleCoords: trend and plunge.
         *      trend  = azimuth of sigma 3 in interval [0, 360), measured clockwise from the North
         *      plunge = vertical angle between the horizontal plane and the sigma 3 axis (positive downward), in interval [-90,90]
         *
         * @param param0
         */
        this.poleS3 = new _types__WEBPACK_IMPORTED_MODULE_0__.PoleCoords();
        // The user selects a Master Principal Stress (MPS) and a Subordinate/Slave Principal Stress (SPS), Sigma 1 or Sigma 3:
        //      The trend and plunge of the MPS are defined by the user
        //      The trend of the SPS is defined by the user while the plunge is calculated from the 3 selected angles
        //      The SPS is located in the plane perpendicular to the MPS.
        this.masterStress = MasterStress.Sigma1;
        this.sigma = [0, 0, 0];
        this.Rrot_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
        this.STP_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
        this.ST_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
        this.poleS1.plunge = plungeS1;
        this.poleS3.plunge = plungeS3;
        this.poleS1.trend = trendS1;
        this.poleS3.trend = trendS3;
        this.sigma = [-1, -stressRatio, 0];
        this.changeMasterStress(masterStress);
        // Perfom the computations...
        this.masterSlave();
        this.rotationTensors(new _types__WEBPACK_IMPORTED_MODULE_0__.SphericalCoords({ phi: this.phiS1, theta: this.thetaS1 }), new _types__WEBPACK_IMPORTED_MODULE_0__.SphericalCoords({ phi: this.phiS3, theta: this.thetaS3 }));
    }
    changeMasterStress(masterStress) {
        this.masterStress = masterStress;
    }
    get ST() {
        return this.ST_;
    }
    get STP() {
        return this.STP_;
    }
    get stressRatio() {
        return this.sigma[1];
    }
    set plungeS1(v) {
        this.poleS1.plunge = v;
        this.masterSlave();
    }
    set trendS1(v) {
        this.poleS1.trend = v;
        this.masterSlave();
    }
    set plungeS3(v) {
        this.poleS3.plunge = v;
        this.masterSlave();
    }
    set trendS3(v) {
        this.poleS3.trend = v;
        this.masterSlave();
    }
    get plunge() {
        if (this.masterStress === MasterStress.Sigma1) {
            return this.poleS3.plunge;
        }
        else {
            return this.poleS1.plunge;
        }
    }
    /**
     * Example:
     * ```ts
     * const s = new StressTensor({trendS1, plungeS1, trendS3, plungeS3, masterStress, sigma})
     * const r = s.Rrot
     * ```
     */
    get Rrot() {
        return this.Rrot_;
    }
    get RTrot() {
        return (0,_types__WEBPACK_IMPORTED_MODULE_0__.transposeTensor)(this.Rrot_);
    }
    masterSlave() {
        // Calculate the plunge of the Slave stress axis
        if (this.masterStress === MasterStress.Sigma1) {
            // The trend and plunge of sigma 1 are set by the user
            const sigma = (0,_types__WEBPACK_IMPORTED_MODULE_0__.lineSphericalCoords)({
                trend: this.poleS1.trend,
                plunge: this.poleS1.plunge
            });
            // The trend of sigma 3 is set by the user while the plunge of sigma 3 has to be calculated
            this.phiS3 = (0,_types__WEBPACK_IMPORTED_MODULE_0__.trend2phi)(this.poleS3.trend);
            this.thetaS3 = this.thetaSlave(sigma.phi, sigma.theta, this.phiS3);
            this.poleS3.plunge = (0,_types__WEBPACK_IMPORTED_MODULE_0__.rad2deg)(this.thetaS3 - Math.PI / 2);
        }
        else {
            // The trend and plunge of sigma 3 are set by the user
            const sigma = (0,_types__WEBPACK_IMPORTED_MODULE_0__.lineSphericalCoords)({
                trend: this.poleS3.trend,
                plunge: this.poleS3.plunge
            });
            // The trend of sigma 1 is set by the user while the plunge of sigma 3 has to be calculated
            this.phiS1 = (0,_types__WEBPACK_IMPORTED_MODULE_0__.trend2phi)(this.poleS1.trend);
            this.thetaS1 = this.thetaSlave(sigma.phi, sigma.theta, this.phiS1);
            this.poleS1.plunge = (0,_types__WEBPACK_IMPORTED_MODULE_0__.rad2deg)(this.thetaS1 - Math.PI / 2);
        }
    }
    thetaSlave(phiMaster, thetaMaster, phiSlave) {
        // thetaSlave: colatitude or polar angle of the slave stress axis measured downward from the zenith (upward direction) [0, PI)
        let thetaSlave = 0;
        if (thetaMaster > 0 && thetaMaster < Math.PI / 2) {
            // The master stress axis is in the upper hemisphere
            // phi_NormalPlane = azimuthal angle of the normal plane in inteerval [0, 2 PI):
            //  The plane is located in the upper hemisphere in the anticlockwise direction realtive to phi_NormalPlane   
            let phi_NormalPlane = phiMaster + Math.PI / 2;
            if (phi_NormalPlane >= 2 * Math.PI) {
                phi_NormalPlane -= 2 * Math.PI;
            }
            // phi_Y_NormalPlane = inclination angle of the normal fault plane (0, PI/2 )()
            let phi_Y_NormalPlane = thetaMaster;
            // phi_NP_Slave = azimuthal angle between phi_NormalPlane and phiSlave
            let phi_NP_Slave = phiSlave - phi_NormalPlane;
            // Analytic equation relating angles in spherical coords 
            thetaSlave = Math.atan(1 / (Math.sin(phi_NP_Slave) * Math.tan(phi_Y_NormalPlane)));
            if (thetaSlave < 0) {
                thetaSlave += Math.PI;
            }
        }
        else if (thetaMaster > Math.PI / 2) {
            // The master stress axis is in the lower hemisphere
            let phi_NormalPlane = phiMaster - Math.PI / 2;
            if (phi_NormalPlane < 0) {
                phi_NormalPlane += 2 * Math.PI;
            }
            // phi_Y_NormalPlane = inclination angle of the normal fault plane (0, PI/2 )()
            let phi_Y_NormalPlane = Math.PI - thetaMaster;
            // phi_NP_Slave = azimuthal angle between phi_NormalPlane and phiSlave
            let phi_NP_Slave = phiSlave - phi_NormalPlane;
            // Analytic equation relating angles in spherical coords 
            thetaSlave = Math.atan(1 / (Math.sin(phi_NP_Slave) * Math.tan(phi_Y_NormalPlane)));
            if (thetaSlave < 0) {
                thetaSlave += Math.PI;
            }
        }
        else if (thetaMaster === 0) {
            // The master stress axis is vertical
        }
        else if (thetaMaster === Math.PI / 2) { // angle in degrees can be more precise for this test
            // The master stress axis is horizontal
        }
        return thetaSlave;
    }
    rotationTensors(sigma_1_sphere, sigma_3_sphere) {
        // Calculate the rotation tensors between reference frames S and S', where
        // S' = (X',Y',Z') is the principal stress reference frame, parallel to (sigma_1, sigma_3, sigma_2);
        // S =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
        //
        // Let Rrot be the rotation tensor R between reference systems S and S', such that:
        //      V' = R V,  where V and V' are the same vector defined in reference frames S and S', respectively
        // The lines of matrix R are given by the unit vectors parallel to X1', X2', and X3' defined in reference system S:
        // Sigma_1 axis: Unit vector e1'. The scalar product: e1'.V = V'(1)
        const Rrot = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
        Rrot[0][0] = Math.sin(sigma_1_sphere.theta) * Math.cos(sigma_1_sphere.phi);
        Rrot[0][1] = Math.sin(sigma_1_sphere.theta) * Math.sin(sigma_1_sphere.phi);
        Rrot[0][2] = Math.cos(sigma_1_sphere.theta);
        // Sigma_3 axis: Unit vector e2'. The scalar product: e2'.V = V'(2)
        Rrot[0][0] = Math.sin(sigma_3_sphere.theta) * Math.cos(sigma_3_sphere.phi);
        Rrot[0][1] = Math.sin(sigma_3_sphere.theta) * Math.sin(sigma_3_sphere.phi);
        Rrot[0][2] = Math.cos(sigma_3_sphere.theta);
        // Sigma_2 axis: Unit vector e3'. The scalar product: e3'.V = V'(3)
        // e3' is calculated from the cross product e3' = e1' x e2' :
        Rrot[2][0] = Rrot[0][1] * Rrot[1][2] - Rrot[0][2] * Rrot[1][1];
        Rrot[2][1] = Rrot[0][2] * Rrot[1][0] - Rrot[0][0] * Rrot[1][2];
        Rrot[2][2] = Rrot[0][0] * Rrot[1][1] - Rrot[0][1] * Rrot[1][0];
        // Let RTrot be the rotation tensor R between reference systems S' and S, such that:
        //      V = RTrot V',  where V and V' are the same vector defined in reference frames S and S', respectively
        this.Rrot_ = Rrot;
        this.stressTensorS();
    }
    stressTensorS() {
        // Calculate the stress tensor ST in reference frame S from the stress tensor in reference frame S':
        //      ST = RTrot STP Rrot
        //
        // where
        //
        //      S' = (X',Y',Z') is the principal stress reference frame, parallel to (sigma_1, sigma_3, sigma_2);
        //      S =  (X, Y, Z ) is the geographic reference frame  oriented in (East, North, Up) directions.
        //      STP = Stress tensor in the principal stress reference frame.
        this.STP_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.stressTensorPrincipalAxes)(this.sigma.map(v => -v));
        const T1 = (0,_types__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: (0,_types__WEBPACK_IMPORTED_MODULE_0__.transposeTensor)(this.Rrot_), B: this.STP });
        this.ST_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors)({ A: T1, B: this.Rrot_ });
    }
}
/*
import { Matrix3x3 } from "./math"
import { PoleCoords } from "./PoleCoords"

export class StressTensor {
    stressTensor_: Matrix3x3
    private rotationTensorR_: Matrix3x3
    private rotationTensorRT_: Matrix3x3
    private polarCoordsS1_: PoleCoords
    private polarCoordsS3_: PoleCoords
    private stressRatio_: number


    constructor(
        stressTensor_: Matrix3x3,
        rotationTensorR_: Matrix3x3,
        rotationTensorRT_: Matrix3x3,
        polarCoordsS1_: PoleCoords,
        polarCoordsS3_: PoleCoords,
        stressRatio_: number)
        {
        // this.stressTensor_ =
        }
    
    get stressTensor(): Matrix3x3 {
        return this.stressTensor_
    }
    set stressTensor(value: Matrix3x3) {
        this.stressTensor_ = value
    }
}
*/


/***/ }),

/***/ "./lib/types/index.ts":
/*!****************************!*\
  !*** ./lib/types/index.ts ***!
  \****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "MasterStress": () => (/* reexport safe */ _StressTensor__WEBPACK_IMPORTED_MODULE_4__.MasterStress),
/* harmony export */   "PoleCoords": () => (/* reexport safe */ _PoleCoords__WEBPACK_IMPORTED_MODULE_2__.PoleCoords),
/* harmony export */   "SphericalCoords": () => (/* reexport safe */ _SphericalCoords__WEBPACK_IMPORTED_MODULE_3__.SphericalCoords),
/* harmony export */   "StressTensor": () => (/* reexport safe */ _StressTensor__WEBPACK_IMPORTED_MODULE_4__.StressTensor),
/* harmony export */   "add_Vectors": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.add_Vectors),
/* harmony export */   "angularDifStriations": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.angularDifStriations),
/* harmony export */   "arcCircle": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.arcCircle),
/* harmony export */   "cloneMatrix3x3": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.cloneMatrix3x3),
/* harmony export */   "constant_x_Vector": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.constant_x_Vector),
/* harmony export */   "crossProduct": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.crossProduct),
/* harmony export */   "deg2rad": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.deg2rad),
/* harmony export */   "faultStressComponents": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.faultStressComponents),
/* harmony export */   "lineSphericalCoords": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.lineSphericalCoords),
/* harmony export */   "minRotAngleRotationTensor": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.minRotAngleRotationTensor),
/* harmony export */   "mohrCircleLine": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.mohrCircleLine),
/* harmony export */   "mohrCirclePoint": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.mohrCirclePoint),
/* harmony export */   "multiplyTensors": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.multiplyTensors),
/* harmony export */   "newMatrix3x3": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.newMatrix3x3),
/* harmony export */   "newMatrix3x3Identity": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.newMatrix3x3Identity),
/* harmony export */   "newPoint3D": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.newPoint3D),
/* harmony export */   "newVector3D": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.newVector3D),
/* harmony export */   "normalVector": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.normalVector),
/* harmony export */   "normalizeVector": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.normalizeVector),
/* harmony export */   "normalizedCrossProduct": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.normalizedCrossProduct),
/* harmony export */   "properRotationTensor": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.properRotationTensor),
/* harmony export */   "rad2deg": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.rad2deg),
/* harmony export */   "rotationParamsFromRotTensor": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.rotationParamsFromRotTensor),
/* harmony export */   "scalarProduct": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.scalarProduct),
/* harmony export */   "scalarProductUnitVectors": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.scalarProductUnitVectors),
/* harmony export */   "setValueInUnitInterval": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.setValueInUnitInterval),
/* harmony export */   "signedAngularDifStriations": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.signedAngularDifStriations),
/* harmony export */   "spherical2unitVectorCartesian": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.spherical2unitVectorCartesian),
/* harmony export */   "stressTensorPrincipalAxes": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.stressTensorPrincipalAxes),
/* harmony export */   "tensor_x_Vector": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.tensor_x_Vector),
/* harmony export */   "transposeTensor": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.transposeTensor),
/* harmony export */   "trend2phi": () => (/* reexport safe */ _mechanics__WEBPACK_IMPORTED_MODULE_6__.trend2phi),
/* harmony export */   "unitVectorCartesian2Spherical": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.unitVectorCartesian2Spherical),
/* harmony export */   "vectorMagnitude": () => (/* reexport safe */ _math__WEBPACK_IMPORTED_MODULE_5__.vectorMagnitude)
/* harmony export */ });
/* harmony import */ var _GenericCurve__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./GenericCurve */ "./lib/types/GenericCurve.ts");
/* harmony import */ var _MohrPoint__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./MohrPoint */ "./lib/types/MohrPoint.ts");
/* harmony import */ var _PoleCoords__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./PoleCoords */ "./lib/types/PoleCoords.ts");
/* harmony import */ var _SphericalCoords__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./SphericalCoords */ "./lib/types/SphericalCoords.ts");
/* harmony import */ var _StressTensor__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./StressTensor */ "./lib/types/StressTensor.ts");
/* harmony import */ var _math__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./math */ "./lib/types/math.ts");
/* harmony import */ var _mechanics__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./mechanics */ "./lib/types/mechanics.ts");









/***/ }),

/***/ "./lib/types/math.ts":
/*!***************************!*\
  !*** ./lib/types/math.ts ***!
  \***************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "add_Vectors": () => (/* binding */ add_Vectors),
/* harmony export */   "cloneMatrix3x3": () => (/* binding */ cloneMatrix3x3),
/* harmony export */   "constant_x_Vector": () => (/* binding */ constant_x_Vector),
/* harmony export */   "crossProduct": () => (/* binding */ crossProduct),
/* harmony export */   "deg2rad": () => (/* binding */ deg2rad),
/* harmony export */   "minRotAngleRotationTensor": () => (/* binding */ minRotAngleRotationTensor),
/* harmony export */   "multiplyTensors": () => (/* binding */ multiplyTensors),
/* harmony export */   "newMatrix3x3": () => (/* binding */ newMatrix3x3),
/* harmony export */   "newMatrix3x3Identity": () => (/* binding */ newMatrix3x3Identity),
/* harmony export */   "newPoint3D": () => (/* binding */ newPoint3D),
/* harmony export */   "newVector3D": () => (/* binding */ newVector3D),
/* harmony export */   "normalVector": () => (/* binding */ normalVector),
/* harmony export */   "normalizeVector": () => (/* binding */ normalizeVector),
/* harmony export */   "normalizedCrossProduct": () => (/* binding */ normalizedCrossProduct),
/* harmony export */   "properRotationTensor": () => (/* binding */ properRotationTensor),
/* harmony export */   "rad2deg": () => (/* binding */ rad2deg),
/* harmony export */   "rotationParamsFromRotTensor": () => (/* binding */ rotationParamsFromRotTensor),
/* harmony export */   "scalarProduct": () => (/* binding */ scalarProduct),
/* harmony export */   "scalarProductUnitVectors": () => (/* binding */ scalarProductUnitVectors),
/* harmony export */   "setValueInUnitInterval": () => (/* binding */ setValueInUnitInterval),
/* harmony export */   "spherical2unitVectorCartesian": () => (/* binding */ spherical2unitVectorCartesian),
/* harmony export */   "stressTensorPrincipalAxes": () => (/* binding */ stressTensorPrincipalAxes),
/* harmony export */   "tensor_x_Vector": () => (/* binding */ tensor_x_Vector),
/* harmony export */   "transposeTensor": () => (/* binding */ transposeTensor),
/* harmony export */   "unitVectorCartesian2Spherical": () => (/* binding */ unitVectorCartesian2Spherical),
/* harmony export */   "vectorMagnitude": () => (/* binding */ vectorMagnitude)
/* harmony export */ });
/* harmony import */ var _SphericalCoords__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./SphericalCoords */ "./lib/types/SphericalCoords.ts");

/**
 * @category Math
 */
function newPoint3D() {
    return [0, 0, 0];
}
/**
 * @category Math
 */
function newVector3D() {
    return [0, 0, 0];
}
/**
 * @category Math
 */
function newMatrix3x3() {
    return [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
}
/**
 * @category Math
 */
function cloneMatrix3x3(m) {
    return [[...m[0]], [...m[1]], [...m[2]]];
}
/**
 * @category Math
 */
function newMatrix3x3Identity() {
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
}
/**
 * @category Math
 */
function deg2rad(a) {
    return a * Math.PI / 180;
}
/**
 * @category Math
 */
const rad2deg = (a) => a / Math.PI * 180;
/**
 * @category Math
 */
function vectorMagnitude(vector) {
    // Calculate the magnitude of the vector
    let magVector = Math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2);
    return magVector;
}
/**
 * @category Math
 */
function normalizeVector(vector, norm) {
    if (norm !== undefined) {
        if (norm === 0) {
            throw new Error(`norm is zero`);
        }
        return [vector[0] / norm, vector[1] / norm, vector[2] / norm];
    }
    // Calculate the magnitude of the vector
    let magVector = vectorMagnitude(vector);
    if (magVector === 0) {
        throw new Error(`vector is null and cannot be normalized`);
    }
    return [vector[0] / magVector, vector[1] / magVector, vector[2] / magVector];
}
/**
 * @category Math
 */
function stressTensorPrincipalAxes(sigma) {
    // Calculate the stress tensor ST in the principal stress frame 
    const STP = newMatrix3x3();
    // Stress tensor in the principal stress axis is diagonal
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            if (i === j) {
                STP[i][j] = sigma[i];
            }
            else {
                STP[i][j] = 0;
            }
        }
    }
    return STP;
}
/**
 * @category Math
 */
function tensor_x_Vector({ T, V }) {
    // Pre-multply tensor T by vector V
    const TV = newVector3D();
    TV[0] = T[0][0] * V[0] + T[0][1] * V[1] + T[0][2] * V[2];
    TV[1] = T[1][0] * V[0] + T[1][1] * V[1] + T[1][2] * V[2];
    TV[2] = T[2][0] * V[0] + T[2][1] * V[1] + T[2][2] * V[2];
    return TV;
}
/**
 * @category Math
 */
function constant_x_Vector({ k, V }) {
    // multiply vector V by constant k
    const TV = newVector3D();
    TV[0] = k * V[0];
    TV[1] = k * V[1];
    TV[2] = k * V[2];
    return TV;
}
/**
 * @category Math
 */
function add_Vectors({ U, V }) {
    // multiply vector V by constant k
    const TV = newVector3D();
    TV[0] = U[0] + V[0];
    TV[1] = U[1] + V[1];
    TV[2] = U[2] + V[2];
    return TV;
}
/**
 * @category Math
 */
function scalarProduct({ U, V }) {
    // Pre-multply tensor T by vector V
    const UdotV = U[0] * V[0] + U[1] * V[1] + U[2] * V[2];
    return UdotV;
}
/**
 * @category Math
 */
function scalarProductUnitVectors({ U, V }) {
    // Pre-multply tensor T by vector V
    let UdotV = U[0] * V[0] + U[1] * V[1] + U[2] * V[2];
    // The scalar product of unit vectors: -1 <= UdotV <= 1
    UdotV = setValueInUnitInterval(UdotV);
    return UdotV;
}
/**
 * @category Math
 * Scalar value is constrained to interval [-1,1]
 */
function setValueInUnitInterval(U) {
    let V = U;
    if (U > 1) {
        V = 1;
    }
    if (U < -1) {
        V = -1;
    }
    return V;
}
/**
 * @brief Calculate the cross product of 2 vectors U and V: U x V
 * @param {U: Vector3, V: Vector3}
 * @example
 * ```ts
 * const Ua: Vector3 = ...
 * const Va: Vector3 = ...
 * const return = crossProduct({U: Ua, V: Va})
 * ```
 * @category Math
 */
function crossProduct({ U, V }) {
    return [U[1] * V[2] - U[2] * V[1],
        U[2] * V[0] - U[0] * V[2],
        U[0] * V[1] - U[1] * V[0]];
}
/**
 * @brief Calculate the cross product of 2 vectors U and V: U x V
 * @param {U: Vector3, V: Vector3}
 * @example
 * ```ts
 * const Ua: Vector3 = ...
 * const Va: Vector3 = ...
 * const return = crossProduct({U: Ua, V: Va})
 * ```
 * @category Math
 */
function normalizedCrossProduct({ U, V }) {
    let W;
    W = crossProduct({ U, V });
    return normalizeVector(W);
}
/**
 * @brief Calculate the multiplication of 2 tensors A and B: Cik = Aij Bjk
 * @param {A: Matrix3x3, B: Matrix3x3}
 * @example
 * ```ts
 * const Ma: Matrix3x3 = ...
 * const Mb: Matrix3x3 = ...
 * const C = multiplyTensors({A: Ma, B: Mb})
 * ```
 * @category Math
 */
function multiplyTensors({ A, B }) {
    const C = newMatrix3x3();
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            C[i][j] = 0;
            for (let k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
/**
 * @category Math
 */
function transposeTensor(A) {
    // Calculate the multiplication of 2 tensors: Cik = Aij Bjk
    const B = newMatrix3x3();
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}
/**
 * @category Math
 */
function rotationParamsFromRotTensor(rotTensor) {
    // The cartesian and spherical coords of a unit vector corresponding to the rotation axis are determined 
    // from the components of the tensor definning a proper rotation
    let rotVector;
    // The axis of rotation is determined form the compoientns of the matrix of a proper rotation
    rotVector[0] = rotTensor[2][1] - rotTensor[1][2];
    rotVector[1] = rotTensor[0][2] - rotTensor[2][0];
    rotVector[2] = rotTensor[1][0] - rotTensor[0][1];
    let rotVectorMag = vectorMagnitude(rotVector);
    // The magnitude of rotVector computed this way is ||rotVector|| = 2 sin θ, where θ is the angle of rotation.
    let rotMag = Math.asin(rotVectorMag / 2);
    let rotAxis = undefined;
    if (rotMag > 0) {
        rotAxis = normalizeVector(rotVector, rotMag);
    }
    else {
        rotAxis = [1, 0, 0];
    }
    return {
        rotAxis,
        rotMag
    };
}
/**
 * @category Math
 */
function normalVector({ phi, theta }) {
    /**
     * Define unit vector normal to the fault plane in the upper hemisphere (pointing upward) from angles in spherical coordinates.
     * The normal vector is constnat for each fault plane and is defined in the geographic reference system: S = (X,Y,Z)
    */
    let normal = newVector3D(); // ***
    normal[0] = Math.sin(phi) * Math.cos(theta);
    normal[1] = Math.sin(phi) * Math.sin(theta);
    normal[2] = Math.cos(phi);
    return normal;
}
/**
 * @category Math
 */
function spherical2unitVectorCartesian(spheriCoords) {
    // The principal stress axes and microfault data such as stylolites can be represented by lines.
    // A line is defined by its trend and plunge angles in the geographic reference system:
    // trend = azimuth of the line in interval [0, 360), measured clockwise from the North direction
    // plunge =  vertical angle between the horizontal plane and the sigma 1 axis (positive downward), in interval [0,90]
    // (phi,theta) : spherical coordinate angles defining the unit vector in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
    // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
    // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
    const V = newVector3D();
    V[0] = Math.sin(spheriCoords.theta) * Math.cos(spheriCoords.phi);
    V[1] = Math.sin(spheriCoords.theta) * Math.sin(spheriCoords.phi);
    V[2] = Math.cos(spheriCoords.theta);
    return V;
}
/**
 * @category Math
 */
function unitVectorCartesian2Spherical(V, EPS = 1e-7) {
    // This routine inverts the following equations in spherical coordinates:
    // V[0] = Math.sin(spheriCoords.theta) * Math.cos(spheriCoords.phi)
    // V[1] = Math.sin(spheriCoords.theta) * Math.sin(spheriCoords.phi)
    // V[2] = Math.cos(spheriCoords.theta)
    // (phi,theta) : spherical coordinate angles defining the unit vector in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
    // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
    // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
    let spheriCoords = new _SphericalCoords__WEBPACK_IMPORTED_MODULE_0__.SphericalCoords();
    // Unit vector component V[0] is constrained to interval [-1,1]
    V[2] = setValueInUnitInterval(V[2]);
    // theta = polar angle in interval [0,PI]
    spheriCoords.theta = Math.acos(V[2]);
    let stheta = Math.sin(spheriCoords.theta);
    if (Math.abs(stheta) > EPS) { // In principle, stheta >=0
        // cphi = cos(phi)
        let cphi = V[0] / stheta;
        // cphi is constrained to interval [-1,1]
        cphi = setValueInUnitInterval(cphi);
        // phi is in interval [0,PI]
        spheriCoords.phi = Math.acos(cphi);
        if (V[1] < 0) {
            // phi is in interval (PI,2PI). The angle is obtained by reflexion on the x axis:
            spheriCoords.phi = 2 * Math.PI - spheriCoords.phi;
        }
    }
    else {
        // theta is close to 0 or PI, thus the unit vector is close to the vertical axis
        // phi can take any value
        spheriCoords.phi = 0;
    }
    return spheriCoords;
}
/**
 * @category Math
 */
function properRotationTensor({ nRot, angle }) {
    // Calculate the proper rotation tensor psi corresponding to a rotation angle around a unit axis nRot
    // Psi allows to calculate the new coords of a vector undergoing a given rotation
    const PsiRot = newMatrix3x3();
    let cosa = Math.cos(angle);
    let sina = Math.sin(angle);
    PsiRot[0][0] = cosa + nRot[0] ** 2 * (1 - cosa);
    PsiRot[0][1] = nRot[0] * nRot[1] * (1 - cosa) - nRot[2] * sina;
    PsiRot[0][2] = nRot[0] * nRot[2] * (1 - cosa) + nRot[1] * sina;
    PsiRot[1][0] = nRot[1] * nRot[0] * (1 - cosa) + nRot[2] * sina;
    PsiRot[1][1] = cosa + nRot[1] ** 2 * (1 - cosa);
    PsiRot[1][2] = nRot[1] * nRot[2] * (1 - cosa) - nRot[0] * sina;
    PsiRot[2][0] = nRot[2] * nRot[0] * (1 - cosa) - nRot[1] * sina;
    PsiRot[2][1] = nRot[2] * nRot[1] * (1 - cosa) + nRot[0] * sina;
    PsiRot[2][2] = cosa + nRot[2] ** 2 * (1 - cosa);
    return PsiRot;
}
/**
 * @category Math
 */
function minRotAngleRotationTensor(rotTensor, EPS = 1e-7) {
    // let rotTensor be the rotation tensor between two right-handed references systems Sa = (Xa, Ya, Za) and Sb = (Xb, Yb, Zb) such that:
    //      Vb = rotTensor Va
    // where Va and Vb are corresponding vectors defined in reference systems Sa and Sb, respectively
    // Sa may correspond to the reference system of the principal directions of a stress tensor defined by a pair of conjugate faults
    // or neoformed striated plane.
    // Sb is the reference system of the principal directions of a stress tensor defined by the interactive search
    // or an inverse method. 
    // Recall that the pricipal stress direction (Sigma 1, Sigma 3, Sigma 2) are parallel to (X, Y, Z), respectively
    // This function calculates the minimum rotation angle between Sa and Sb, by considering the four possible right-handed reference systems
    // that are consistent with pricipal stress directions en system Sb. 
    // The angle of rotation associated to rotTensor is defined by the trace tr(rotTensor), according to the relation:
    //      tr(rotTensor) = 1 + 2 cos(rotAngle)
    // where rotAngle is in interval [0,PI]
    // Note that the inverse rotation tensor defined by the transpoded matrix has the same trace. 
    // Thus the rotation angle is the same for tensors rotTensor and rotTensorT (i.e., transposed)
    let traceRotTensor = new Array(4);
    // The trace of the first rotation tensor such that reference system Sb0 = (Xb, Yb, Zb)
    traceRotTensor[0] = rotTensor[0][0] + rotTensor[1][1] + rotTensor[2][2];
    // The trace of the second rotation tensor such that reference system Sb1 = (Xb, -Yb, -Zb)
    // System Sb1 is obtained by rotating Sb0 at an angle of PI around Xb
    // Note that Sb1 is right-handed and its principal axes are parallel to (Sigma 1, Sigma 3, Sigma 2)
    traceRotTensor[1] = rotTensor[0][0] - rotTensor[1][1] - rotTensor[2][2];
    // The trace of the second rotation tensor such that reference system Sb2 = (-Xb, Yb, -Zb)
    // System Sb2 is obtained by rotating Sb0 at an angle of PI around Yb
    // Note that Sb2 is right-handed and its principal axes are parallel to (Sigma 1, Sigma 3, Sigma 2)
    traceRotTensor[2] = -rotTensor[0][0] + rotTensor[1][1] - rotTensor[2][2];
    // The trace of the second rotation tensor such that reference system Sb3 = (Xb, -Yb, -Zb)
    // System Sb3 is obtained by rotating Sb0 at an angle of PI around Zb
    // Note that Sb3 is right-handed and its principal axes are parallel to (Sigma 1, Sigma 3, Sigma 2)
    traceRotTensor[3] = -rotTensor[0][0] - rotTensor[1][1] + rotTensor[2][2];
    const max = traceRotTensor.reduce((cur, v) => Math.max(cur, v), Number.NEGATIVE_INFINITY);
    let cosMinRotAngle = (max - 1) / 2;
    if (Math.abs(cosMinRotAngle) > 1) {
        if (Math.abs(cosMinRotAngle) > 1 + EPS) {
            throw new Error(`The cosine of the minimum rotation angle of the rotation tensor is not in the unit interval`);
        }
        cosMinRotAngle = setValueInUnitInterval(cosMinRotAngle);
    }
    return Math.acos(cosMinRotAngle);
}


/***/ }),

/***/ "./lib/types/mechanics.ts":
/*!********************************!*\
  !*** ./lib/types/mechanics.ts ***!
  \********************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "angularDifStriations": () => (/* binding */ angularDifStriations),
/* harmony export */   "arcCircle": () => (/* binding */ arcCircle),
/* harmony export */   "faultStressComponents": () => (/* binding */ faultStressComponents),
/* harmony export */   "lineSphericalCoords": () => (/* binding */ lineSphericalCoords),
/* harmony export */   "mohrCircleLine": () => (/* binding */ mohrCircleLine),
/* harmony export */   "mohrCirclePoint": () => (/* binding */ mohrCirclePoint),
/* harmony export */   "signedAngularDifStriations": () => (/* binding */ signedAngularDifStriations),
/* harmony export */   "trend2phi": () => (/* binding */ trend2phi)
/* harmony export */ });
/* harmony import */ var _math__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./math */ "./lib/types/math.ts");
/* harmony import */ var _SphericalCoords__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./SphericalCoords */ "./lib/types/SphericalCoords.ts");
/* harmony import */ var _analysis_Curve3D__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../analysis/Curve3D */ "./lib/analysis/Curve3D.ts");



/**
 * @category Mechanics
 */
function faultStressComponents({ stressTensor, normal }) {
    // Calculate the stress components applied on a fault plane as a result of a stress tensor StressTensor defined in the reference system S
    // normal: unit vector normal to the fault plane (pointing upward) defined in the geographic reference system: S = (X,Y,Z)
    // Calculate total stress vector
    let stress = (0,_math__WEBPACK_IMPORTED_MODULE_0__.tensor_x_Vector)({ T: stressTensor, V: normal });
    // Calculate normal stress (positive = extension, negative = compression). 
    // In principle the normal stress is negative since the principal stresses are <= 0.
    let normalStress = (0,_math__WEBPACK_IMPORTED_MODULE_0__.scalarProduct)({ U: stress, V: normal });
    // Calculate the shear stress vector in reference system S
    // shearStress + normalStress * normal = stress (i.e., total stress)
    let shearStress = [0, 0, 0];
    shearStress[0] = stress[0] - normalStress * normal[0];
    shearStress[1] = stress[1] - normalStress * normal[1];
    shearStress[2] = stress[2] - normalStress * normal[2];
    let shearStressMag = (0,_math__WEBPACK_IMPORTED_MODULE_0__.vectorMagnitude)(shearStress);
    return {
        shearStress,
        normalStress,
        shearStressMag
    };
}
/**
 * @category Mechanics
 */
function angularDifStriations({ e_striation, shearStress, shearStressMag }) {
    // Calculate the angular difference between the measured and calculated striations
    // The angular difference calculated by the scalar product is unsigned (in interval [0,Pi])
    let angularDifStriae = 0;
    if (shearStressMag > 0) {
        // The angular difference is calculated using the scalar product: 
        //      e_striation . shearStress = |e_striation| |shearStress| cos(angularDifStriae) = 1 . shearStressMag . cos(angularDifStriae)
        angularDifStriae = (Math.acos(e_striation[0] * shearStress[0] + e_striation[1] * shearStress[1]
            + e_striation[2] * shearStress[2])) / shearStressMag;
    }
    else {
        // The calculated shear stress is zero (i.e., the fault plane is parallel to a principal stress)
        // In such situation we may consider that the calculated striation can have any direction.
    }
    return angularDifStriae;
}
/**
 * @category Mechanics
 */
function signedAngularDifStriations({ normal, e_striation, shearStress, shearStressMag }) {
    // Calculate the angular difference between the measured and calculated striations
    // The angular difference calculated by the scalar product is unsigned (in interval [0,Pi])
    let angularDifStriae = 0;
    let nAux;
    if (shearStressMag > 0) {
        // The angular difference is calculated using the scalar product in interval [0,PI]: 
        //      e_striation . shearStress = |e_striation| |shearStress| cos(angularDifStriae) = 1 . shearStressMag . cos(angularDifStriae)
        angularDifStriae = (Math.acos(e_striation[0] * shearStress[0] + e_striation[1] * shearStress[1]
            + e_striation[2] * shearStress[2])) / shearStressMag;
        // nAux = auxiliary normal vector pointing outward (i.e., in the same direction as normal vector) if the angle between the measured and calculated striation
        //        is positive - anti-clockwise in interval (0,PI).
        //        Conversely, nAux points in the opposite direction if the angle is clockwise - negative, in interval (-PI,0)           
        let nAux = (0,_math__WEBPACK_IMPORTED_MODULE_0__.crossProduct)({ U: e_striation, V: shearStress });
        let sAux = (0,_math__WEBPACK_IMPORTED_MODULE_0__.scalarProduct)({ U: normal, V: nAux });
        if (sAux < 0) {
            // angularDifStriae is negative (i.e., clokwise, in interval (-PI,0) )
            angularDifStriae = -angularDifStriae;
        }
    }
    else {
        // The calculated shear stress is zero (i.e., the fault plane is parallel to a principal stress)
        // In such situation we may consider that the calculated striation can have any direction.
    }
    return angularDifStriae;
}
/**
 * @category Mechanics
 */
function mohrCircleLine({ r, first, second, sigma_1, sigma_2, sigma_3 }) {
    const lineBuilder = new _analysis_Curve3D__WEBPACK_IMPORTED_MODULE_2__.Curve3D();
    if (sigma_2 == sigma_3) {
        // Particular Case 1: revolution stress tensor around sigma_1
        // In such situation, points pO and p1 have equal coordinates and alfa angles
        // Curve is a circle sweeping at an angle alfa0 around sigma_1
        return arcCircle({ r: r, sigma: 'sigma_1', alpha: first.angle });
    }
    else if (sigma_2 == sigma_1) {
        // Particular Case 2: revolution stress tensor around sigma_3
        // In such situation, points pO and p1 have equal coordinates and alfa angles
        // Curve is a circle sweeping at an angle PI/2 - alfa0 around sigma_3
        return arcCircle({ r: r, sigma: 'sigma_3', alpha: Math.PI / 2 - first.angle });
    }
    else {
        // General Case:
        // Add to the list the initial point of the line segment located in one of the 3 Mohr circles
        lineBuilder.addPoint(mohrCirclePoint({ mohrCirc: first.circle, alfa: first.angle, r }));
        // Add to the list the intermediate points of the line segment located between 2 Mohr circles
        // We calculate the direction cosines of the unit vector normal to the fault whose stress state is given by (X, Y)
        //      Note that (X,Y) are the normal and shear stress of a moving point along the Mohr-Coulomb line segment 
        // Without loss of generality, we suppose a stress tensor in strike-slip regime (fixing sigma_1 Eastward and sigma_3 Northward)
        let sigma_X = sigma_1;
        let sigma_Y = sigma_3;
        let sigma_Z = sigma_2;
        for (let i = 1; i < 180; ++i) {
            let X = first.p[0] + (second.p[0] - first.p[0]) * i / 180;
            let Y = first.p[1] + (second.p[1] - first.p[1]) * i / 180;
            let nx = Math.sqrt(((sigma_Y - X) * (sigma_Z - X) + Y ** 2) / ((sigma_Y - sigma_X) * (sigma_Z - sigma_X)));
            let ny = Math.sqrt(((sigma_Z - X) * (sigma_X - X) + Y ** 2) / ((sigma_Z - sigma_Y) * (sigma_X - sigma_Y)));
            let nz = Math.sqrt(((sigma_X - X) * (sigma_Y - X) + Y ** 2) / ((sigma_X - sigma_Z) * (sigma_Y - sigma_Z)));
            const x = r * nx;
            const y = r * ny;
            const z = r * nz;
            lineBuilder.addPoint(x, y, z);
        }
        // Add to the list the final point of the line segment located in one of the 3 Mohr circles
        lineBuilder.addPoint(mohrCirclePoint({ mohrCirc: second.circle, alfa: second.angle, r }));
        return lineBuilder.buffer;
    }
}
/**
 * @category Mechanics
 */
function mohrCirclePoint({ r, mohrCirc, alfa }) {
    // Add to the list the initial or final point of the line segment located in one of the 3 Mohr circles
    if (mohrCirc == '3_1') {
        // The point is located in Mohr circle between sigma_3 and sigma_1
        // alfa is the azimuthal angle in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
        const theta = Math.PI / 2;
        const phi = alfa;
        return [
            r * Math.sin(theta) * Math.cos(phi),
            r * Math.sin(theta) * Math.sin(phi),
            r * Math.cos(theta)
        ];
    }
    else if (mohrCirc == '3_2') {
        // The point is located in Mohr circle between sigma_3 and sigma_2
        // alfa is the polar angle in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
        const theta = alfa;
        const phi = Math.PI / 2;
        return [
            r * Math.sin(theta) * Math.cos(phi),
            r * Math.sin(theta) * Math.sin(phi),
            r * Math.cos(theta)
        ];
    }
    else {
        // MohrCirc == '2_1'
        // The point is located in Mohr circle between sigma_2 and sigma_1
        // alfa is the latitude angle in the reference frame (x,y,z) = (sigma_1,sigma_3, sigma_2) = (East, North, Up)
        const theta = Math.PI / 2 - alfa;
        const phi = 0;
        return [
            r * Math.sin(theta) * Math.cos(phi),
            r * Math.sin(theta) * Math.sin(phi),
            r * Math.cos(theta)
        ];
    }
}
/**
 * Usage:
 * ```ts
 * const arc = arcCircle({alpha: deg2rad(12), sigma: '3_1'})
 * ```
 * @category Mechanics
 */
function arcCircle({ r, sigma, alpha }) {
    const lineBuilder = new _analysis_Curve3D__WEBPACK_IMPORTED_MODULE_2__.Curve3D();
    if (sigma === 'sigma_1') {
        // Particular case 1: sigma_2 = sigma_3 (revolution stress tensor around sigma_1)
        // Generate a circular segment (one quarter of a circle) sweeping an angle alfa around sigma_1
        const x = r * Math.cos(alpha);
        const rad_circle = r * Math.sin(alpha);
        for (let i = 1; i <= 180; ++i) {
            let beta = Math.PI * i / 360;
            const y = rad_circle * Math.cos(beta);
            const z = rad_circle * Math.sin(beta);
            lineBuilder.addPoint(x, y, z);
        }
    }
    else if (sigma === 'sigma_3') {
        // Particular case 2: sigma_2 = sigma_1 (revolution stress tensor around sigma_3)
        // Generate a circular segment (one quarter of a circle) sweeping an angle alfa around sigma_3
        const y = r * Math.cos(alpha);
        const rad_circle = r * Math.sin(alpha);
        for (let i = 1; i <= 180; ++i) {
            let beta = Math.PI * i / 360;
            const x = rad_circle * Math.cos(beta);
            const z = rad_circle * Math.sin(beta);
            lineBuilder.addPoint(x, y, z);
        }
    }
    return lineBuilder.buffer;
}
/**
 * @category Mechanics
 */
function lineSphericalCoords({ trend, plunge }) {
    // The principal stress axes and microfault data such as stylolites can be represented by lines.
    // A line is defined by its trend and plunge angles in the geographic reference system:
    // trend = azimuth of the line in interval [0, 360), measured clockwise from the North direction
    // plunge =  vertical angle between the horizontal plane and the sigma 1 axis (positive downward), in interval [0,90]
    // (phi,theta) : spherical coordinate angles defining the unit vector in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
    // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
    // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
    let phi = trend2phi(trend);
    let theta = (0,_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(plunge) + Math.PI / 2;
    return new _SphericalCoords__WEBPACK_IMPORTED_MODULE_1__.SphericalCoords({ phi, theta });
}
/**
 * @category Mechanics
 */
function trend2phi(trend) {
    // Calculate the value of phi from the trend
    // trend = azimuth of the line in interval [0, 360), measured clockwise from the North direction
    // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
    // trend + phi = 5 PI / 2
    let phi = 5 * Math.PI / 2 - (0,_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(trend);
    if (phi >= 2 * Math.PI) {
        phi -= 2 * Math.PI;
    }
    return phi;
}


/***/ }),

/***/ "./lib/utils/CompactionShearBands.ts":
/*!*******************************************!*\
  !*** ./lib/utils/CompactionShearBands.ts ***!
  \*******************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "CompactionShearBands": () => (/* binding */ CompactionShearBands)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _Fault__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Fault */ "./lib/utils/Fault.ts");


class CompactionShearBands {
    constructor() {
        // (e_phi, e_theta) = unit vectors defining local reference frame tangent to the sphere in spherical coordinates
        this.e_phi = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.e_theta = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.normal_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.nStriation = undefined;
        this.EPS = 1e-7;
    }
    static create({ strike, dipDirection, dip, sensOfMovement, rake, strikeDirection }) {
        const f = new _Fault__WEBPACK_IMPORTED_MODULE_1__.Fault({ strike, dipDirection, dip });
        f.setStriation({ sensOfMovement, rake, strikeDirection });
        return {
            nPlane: f.normal,
            nStriation: f.striation,
            nPerpStriation: f.e_perp_striation,
            fault: f
        };
    }
    get normal() {
        return this.normal_;
    }
    compactionShearBandCheckMouvement({ noPlane, nPlane, coordinates, sensOfMovement, nSigma3_Sc, nSigma2_Sc }) {
        // Function calculating the striation unit vector in the local reference frame in polar coordinates from the rake
        //      The effect of fault movement on the striation is considered in function faultStriationAngle_B
        // Each fault is defined by a set of parameters as follows:
        //      The fault plane orientation is defined by three parameters:
        //      Fault strike: clockwise angle measured from the North direction [0, 360)
        //      Strike direction (optional): (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Fault dip: [0, 90]
        //      Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
        // The orientation of the striation in the fault plane can defined in two different ways (which are exclusive):
        // 1-   Rake (or pitch) [0,90], measured from the strike direction, which points in one of the two opposite directions of the fault strike.
        //      Strike direction : (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Note that the specified strike direction is used to determine the spatial orientation of the striation 
        // 2-   For shallow-dipping planes (i.e., the compass inclinometer is inaccurate):
        //      Striae trend: [0, 360)
        // alphaStria : striation angle measured in the local reference plane (e_phi, e_theta) indicating the motion of the outward block
        //      alphaStria is measured clockwise from e_phi, in interval [0, 2 PI) (this choice is consistent with the definition of the rake, which is measured from the fault strike)
        // (coordinates.phi, coordinates.theta) : spherical coords of compaction-shear band in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
        // Sperical coords are calculated using method faultSphericalCoords in class fault
        // e_phi : unit vector pointing toward the azimuthal direction in the local reference frame in spherical coords
        this.e_phi[0] = -Math.sin(coordinates.phi);
        this.e_phi[1] = Math.cos(coordinates.phi);
        this.e_phi[2] = 0;
        // e_theta : unit vector pointing toward the dip direction in the local reference frame in spherical coords
        this.e_theta[0] = Math.cos(coordinates.theta) * Math.cos(coordinates.phi);
        this.e_theta[1] = Math.cos(coordinates.theta) * Math.sin(coordinates.phi);
        this.e_theta[2] = -Math.sin(coordinates.theta);
        // Check that the sense of mouvement is consistent with the orientation of stress axes
        // This requires the calculation of the striation vector indicating movement of the outward block relative to the inner block
        // The striation vector is in the plane of movement:  nStriation = nPlane x nSigma2_Sc
        let nStriation = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: nPlane, V: nSigma2_Sc });
        if ((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nPlane, V: nSigma3_Sc }) < 0) {
            // nSigma3_Sc points inward (toward the fault plane)
            // Thus, invert the sense of nSigma3_Sc so that it points outward (in the direction of the normal to the plane):
            nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -1, V: nSigma3_Sc });
        }
        if ((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: nSigma3_Sc }) < 0) {
            // nSigma3_Sc and nStriation should be both located in the compressional quadrant relative to the outward hemisphere of the fault plane
            //      In other words, the angle (nSigma3_Sc, nStriation) < PI/2
            // However, if nStriation . nSigma3_Sc < 0 then this condition is not satisfied (i.e. nStriation is located in the dilatant quadrant)
            // Thus, invert the sense of the striation so that it points toward the compressional quadrant:
            nStriation = (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -1, V: nStriation });
        }
        // The strike slip component of movement is defined by the projection of the striation vector along unit vector e_phi
        let strikeSlipComponent = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: this.e_phi });
        // The dip component of movement is defined by the projection of the striation vector along unit vector e_theta
        let dipComponent = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: this.e_theta });
        // Check consistency of strike-lateral component of movement
        if (strikeSlipComponent > this.EPS) {
            // In principle, the conjugated plane has a left-lateral component of movement
            if (sensOfMovement === 3 /* SensOfMovement.RL */ || sensOfMovement === 5 /* SensOfMovement.N_RL */ || sensOfMovement === 7 /* SensOfMovement.I_RL */) {
                // throw new Error('Sense of movement of conjugated fault ' + noPlane + ' includes a right-lateral (RL) which is not consistent with fault kenmatics')
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a right-lateral (RL) which is not consistent with fault kenmatics`);
            }
        }
        else if (strikeSlipComponent < -this.EPS) {
            // In principle, the conjugated plane has a right-lateral component of movement
            if (sensOfMovement === 4 /* SensOfMovement.LL */ || sensOfMovement === 6 /* SensOfMovement.N_LL */ || sensOfMovement === 8 /* SensOfMovement.I_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a left-lateral (LL), which is not consistent with fault kenmatics`);
            }
        }
        else {
            // In principle, the strike-slip component of movement of the conjugated plane is negligeable
            if (sensOfMovement !== 1 /* SensOfMovement.N */ && sensOfMovement !== 2 /* SensOfMovement.I */ && sensOfMovement !== 9 /* SensOfMovement.UKN */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a strike-slip component, which is not consistent with fault kenmatics`);
            }
        }
        // Check consistency of dip-slip component of movement
        if (dipComponent > this.EPS) {
            // In principle, the conjugated plane has a normal component of movement
            if (sensOfMovement === 2 /* SensOfMovement.I */ || sensOfMovement === 7 /* SensOfMovement.I_RL */ || sensOfMovement === 8 /* SensOfMovement.I_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes an inverse (I) component, which is not consistent with fault kenmatics`);
            }
        }
        else if (dipComponent < -this.EPS) {
            // In principle, the conjugated plane has an inverse component of movement
            if (sensOfMovement === 1 /* SensOfMovement.N */ || sensOfMovement === 5 /* SensOfMovement.N_RL */ || sensOfMovement === 6 /* SensOfMovement.N_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a normal component, which is not consistent with fault kenmatics`);
            }
        }
        else {
            // In principle, the dip component of movement of the conjugated plane is negligeable
            if (sensOfMovement !== 3 /* SensOfMovement.RL */ && sensOfMovement !== 4 /* SensOfMovement.LL */ && sensOfMovement !== 9 /* SensOfMovement.UKN */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes an dip component (I or N), which is not consistent with fault kenmatics`);
            }
        }
    }
}


/***/ }),

/***/ "./lib/utils/ConjugateFaults.ts":
/*!**************************************!*\
  !*** ./lib/utils/ConjugateFaults.ts ***!
  \**************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "ConjugateFaults": () => (/* binding */ ConjugateFaults)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");
/* harmony import */ var _Fault__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Fault */ "./lib/utils/Fault.ts");


class ConjugateFaults {
    constructor() {
        // (e_phi, e_theta) = unit vectors defining local reference frame tangent to the sphere in spherical coordinates
        this.e_phi = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.e_theta = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.normal_ = (0,_types__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.nStriation = undefined;
        this.EPS = 1e-7;
    }
    static create({ strike, dipDirection, dip, sensOfMovement, rake, strikeDirection }) {
        const f = new _Fault__WEBPACK_IMPORTED_MODULE_1__.Fault({ strike, dipDirection, dip });
        f.setStriation({ sensOfMovement, rake, strikeDirection });
        return {
            nPlane: f.normal,
            nStriation: f.striation,
            nPerpStriation: f.e_perp_striation,
            fault: f
        };
    }
    get normal() {
        return this.normal_;
    }
    conjugatePlaneCheckMouvement({ noPlane, nPlane, coordinates, sensOfMovement, nSigma3_Sc, nSigma2_Sc }) {
        // Function calculating the striation unit vector in the local reference frame in polar coordinates from the rake
        //      The effect of fault movement on the striation is considered in function faultStriationAngle_B
        // Each fault is defined by a set of parameters as follows:
        //      The fault plane orientation is defined by three parameters:
        //      Fault strike: clockwise angle measured from the North direction [0, 360)
        //      Strike direction (optional): (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Fault dip: [0, 90]
        //      Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
        // The orientation of the striation in the fault plane can defined in two different ways (which are exclusive):
        // 1-   Rake (or pitch) [0,90], measured from the strike direction, which points in one of the two opposite directions of the fault strike.
        //      Strike direction : (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Note that the specified strike direction is used to determine the spatial orientation of the striation 
        // 2-   For shallow-dipping planes (i.e., the compass inclinometer is inaccurate):
        //      Striae trend: [0, 360)
        // alphaStria : striation angle measured in the local reference plane (e_phi, e_theta) indicating the motion of the outward block
        //      alphaStria is measured clockwise from e_phi, in interval [0, 2 PI) (this choice is consistent with the definition of the rake, which is measured from the fault strike)
        // (coordinates.phi, coordinates.theta) : spherical coords of conjugate faults in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
        // Sperical coords are calculated using method faultSphericalCoords in class fault
        // e_phi : unit vector pointing toward the azimuthal direction in the local reference frame in spherical coords
        this.e_phi[0] = -Math.sin(coordinates.phi);
        this.e_phi[1] = Math.cos(coordinates.phi);
        this.e_phi[2] = 0;
        // e_theta : unit vector pointing toward the dip direction in the local reference frame in spherical coords
        this.e_theta[0] = Math.cos(coordinates.theta) * Math.cos(coordinates.phi);
        this.e_theta[1] = Math.cos(coordinates.theta) * Math.sin(coordinates.phi);
        this.e_theta[2] = -Math.sin(coordinates.theta);
        // Check that the sense of mouvement is consistent with the orientation of stress axes
        // This requires the calculation of the striation vector indicating movement of the outward block relative to the inner block
        // The striation vector is in the plane of movement:  nStriation = nPlane x nSigma2_Sc
        let nStriation = (0,_types__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct)({ U: nPlane, V: nSigma2_Sc });
        if ((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nPlane, V: nSigma3_Sc }) < 0) {
            // nSigma3_Sc points inward (toward the fault plane)
            // Thus, invert the sense of nSigma3_Sc so that it points outward (in the direction of the normal to the plane):
            nSigma3_Sc = (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -1, V: nSigma3_Sc });
        }
        if ((0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: nSigma3_Sc }) < 0) {
            // nSigma3_Sc and nStriation should be both located in the compressional quadrant relative to the outward hemisphere of the fault plane
            //      In other words, the angle (nSigma3_Sc, nStriation) < PI/2
            // However, if nStriation . nSigma3_Sc < 0 then this condition is not satisfied (i.e. nStriation is located in the dilatant quadrant)
            // Thus, invert the sense of the striation so that it points toward the compressional quadrant:
            nStriation = (0,_types__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector)({ k: -1, V: nStriation });
        }
        // The strike-slip component of movement is defined by the projection of the striation vector along unit vector e_phi
        let strikeSlipComponent = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: this.e_phi });
        // The dip component of movement is defined by the projection of the striation vector along unit vector e_theta
        let dipComponent = (0,_types__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors)({ U: nStriation, V: this.e_theta });
        // Check consistency of strike-lateral component of movement
        if (strikeSlipComponent > this.EPS) {
            // In principle, the conjugated plane has a left-lateral component of movement
            if (sensOfMovement === 3 /* SensOfMovement.RL */ || sensOfMovement === 5 /* SensOfMovement.N_RL */ || sensOfMovement === 7 /* SensOfMovement.I_RL */) {
                // throw new Error('Sense of movement of conjugated fault ' + noPlane + ' includes a right-lateral (RL) which is not consistent with fault kenmatics')
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a right-lateral (RL) which is not consistent with fault kenmatics`);
            }
        }
        else if (strikeSlipComponent < -this.EPS) {
            // In principle, the conjugated plane has a right-lateral component of movement
            if (sensOfMovement === 4 /* SensOfMovement.LL */ || sensOfMovement === 6 /* SensOfMovement.N_LL */ || sensOfMovement === 8 /* SensOfMovement.I_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a left-lateral (LL), which is not consistent with fault kenmatics`);
            }
        }
        else {
            // In principle, the strike-slip component of movement of the conjugated plane is negligeable
            if (sensOfMovement !== 1 /* SensOfMovement.N */ && sensOfMovement !== 2 /* SensOfMovement.I */ && sensOfMovement !== 9 /* SensOfMovement.UKN */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a strike-slip component, which is not consistent with fault kenmatics`);
            }
        }
        // Check consistency of dip-slip component of movement
        if (dipComponent > this.EPS) {
            // In principle, the conjugated plane has a normal component of movement
            if (sensOfMovement === 2 /* SensOfMovement.I */ || sensOfMovement === 7 /* SensOfMovement.I_RL */ || sensOfMovement === 8 /* SensOfMovement.I_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes an inverse (I) component, which is not consistent with fault kenmatics`);
            }
        }
        else if (dipComponent < -this.EPS) {
            // In principle, the conjugated plane has an inverse component of movement
            if (sensOfMovement === 1 /* SensOfMovement.N */ || sensOfMovement === 5 /* SensOfMovement.N_RL */ || sensOfMovement === 6 /* SensOfMovement.N_LL */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes a normal component, which is not consistent with fault kenmatics`);
            }
        }
        else {
            // In principle, the dip component of movement of the conjugated plane is negligeable
            if (sensOfMovement !== 3 /* SensOfMovement.RL */ && sensOfMovement !== 4 /* SensOfMovement.LL */ && sensOfMovement !== 9 /* SensOfMovement.UKN */) {
                throw new Error(`Sense of movement of conjugated fault ${noPlane} includes an dip component (I or N), which is not consistent with fault kenmatics`);
            }
        }
    }
}


/***/ }),

/***/ "./lib/utils/Fault.ts":
/*!****************************!*\
  !*** ./lib/utils/Fault.ts ***!
  \****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "Direction": () => (/* binding */ Direction),
/* harmony export */   "Fault": () => (/* binding */ Fault),
/* harmony export */   "SensOfMovement": () => (/* binding */ SensOfMovement),
/* harmony export */   "directionExists": () => (/* binding */ directionExists),
/* harmony export */   "faultParams": () => (/* binding */ faultParams),
/* harmony export */   "getDirectionFromString": () => (/* binding */ getDirectionFromString),
/* harmony export */   "getSensOfMovementFromString": () => (/* binding */ getSensOfMovementFromString)
/* harmony export */ });
/* harmony import */ var _types_math__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types/math */ "./lib/types/math.ts");
/* harmony import */ var _types_SphericalCoords__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types/SphericalCoords */ "./lib/types/SphericalCoords.ts");
// Calculate the stress components of fault planes



/**
 * Usage:
 * ```ts
 * const sens = SensOfMovement.LL
 * ```
 * @category Data
 */
var SensOfMovement;
(function (SensOfMovement) {
    SensOfMovement[SensOfMovement["N"] = 1] = "N";
    SensOfMovement[SensOfMovement["I"] = 2] = "I";
    SensOfMovement[SensOfMovement["RL"] = 3] = "RL";
    SensOfMovement[SensOfMovement["LL"] = 4] = "LL";
    SensOfMovement[SensOfMovement["N_RL"] = 5] = "N_RL";
    SensOfMovement[SensOfMovement["N_LL"] = 6] = "N_LL";
    SensOfMovement[SensOfMovement["I_RL"] = 7] = "I_RL";
    SensOfMovement[SensOfMovement["I_LL"] = 8] = "I_LL";
    SensOfMovement[SensOfMovement["UKN"] = 9] = "UKN";
})(SensOfMovement || (SensOfMovement = {}));
function getSensOfMovementFromString(s) {
    switch (s) {
        case 'N': return 1 /* SensOfMovement.N */;
        case 'I': return 2 /* SensOfMovement.I */;
        case 'RL': return 3 /* SensOfMovement.RL */;
        case 'LL': return 4 /* SensOfMovement.LL */;
        case 'N_RL': return 5 /* SensOfMovement.N_RL */;
        case 'N_LL': return 6 /* SensOfMovement.N_LL */;
        case 'I_RL': return 7 /* SensOfMovement.I_RL */;
        case 'I_LL': return 8 /* SensOfMovement.I_LL */;
        case 'UKN': return 9 /* SensOfMovement.UKN */;
    }
}
/**
 * @category Data
 */
var Direction;
(function (Direction) {
    Direction[Direction["E"] = 0] = "E";
    Direction[Direction["W"] = 1] = "W";
    Direction[Direction["N"] = 2] = "N";
    Direction[Direction["S"] = 3] = "S";
    Direction[Direction["NE"] = 4] = "NE";
    Direction[Direction["SE"] = 5] = "SE";
    Direction[Direction["SW"] = 6] = "SW";
    Direction[Direction["NW"] = 7] = "NW";
    Direction[Direction["UNKOWN"] = 8] = "UNKOWN";
    Direction[Direction["ERROR"] = 9] = "ERROR";
})(Direction || (Direction = {}));
const dirs = ['E', 'W', 'N', 'S', 'NE', 'SE', 'SW', 'NW'];
function directionExists(s) {
    if (s.length === 0) {
        return true;
    }
    return dirs.includes(s);
}
function getDirectionFromString(s) {
    if (s.length === 0) {
        return 8 /* Direction.UNKOWN */;
    }
    switch (s) {
        case 'E': return 0 /* Direction.E */;
        case 'W': return 1 /* Direction.W */;
        case 'N': return 2 /* Direction.N */;
        case 'S': return 3 /* Direction.S */;
        case 'NE': return 4 /* Direction.NE */;
        case 'SE': return 5 /* Direction.SE */;
        case 'SW': return 6 /* Direction.SW */;
        case 'NW': return 7 /* Direction.NW */;
        default: return 9 /* Direction.ERROR */;
    }
}
function faultParams({ strike, dipDirection, dip }) {
}
// Each fault comprises information concerning the geometry, the stress parameters and the kinematic parameters:
/**
 * A fault is represented by a plane, this with a normal and a position.
 *
 * Usage:
 * ```ts
 * const f = new Fault({strike: 30, dipDirection: Direction.E, dip: 60})
 * f.setStriation({rake: 20, strikeDirection: Direction.N, sensMouv: 'LL'})
 * ```
 */
class Fault {
    static create({ strike, dipDirection, dip, sensOfMovement, rake, strikeDirection }) {
        const f = new Fault({ strike, dipDirection, dip });
        if (strikeDirection !== undefined) {
            f.setStriation({ sensOfMovement, rake, strikeDirection });
            return {
                nPlane: f.normal,
                nStriation: f.striation,
                nPerpStriation: f.e_perp_striation
            };
        }
        return {
            nPlane: f.normal
        };
    }
    get sphericalCoords() {
        return this.coordinates;
    }
    get normal() {
        return this.normal_;
    }
    /**
     * @brief Get the striation vector in reference system S
     */
    get striation() {
        return this.e_striation_;
    }
    /**
     * @brief Get the striation vector in reference system S
     */
    get e_perp_striation() {
        return this.e_perp_striation_;
    }
    /**
     * Set the orientation of the striation in the fault plane, which can defined in two different ways (which are exclusive):
     * 1. Rake (or pitch) [0,90], measured from the strike direction, which points in one of the two opposite directions of the fault strike.
     *   Strike direction : (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
     * 2. For shallow-dipping planes (i.e., the compass inclinometer is inaccurate):
     *   Striae trend: [0, 360)
     */
    constructor({ strike, dipDirection, dip }) {
        // We define 2 orthonormal right-handed reference systems:
        //      S =  (X, Y, Z ) is the geographic reference frame oriented in (East, North, Up) directions.
        //      S' = (X',Y',Z') is the principal stress reference frame, parallel to (sigma_1, sigma_3, sigma_2);
        // (phi,theta) : spherical coordinate angles defining the unit vector perpendicular to the fault plane (pointing upward)
        //                 in reference system S
        // phi : azimuth phi in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
        // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
        this.coordinates = new _types_SphericalCoords__WEBPACK_IMPORTED_MODULE_1__.SphericalCoords();
        // private phi:    number      // constant value for each fault plane
        // private theta:  number      // constant value for each fault plane
        // normal: unit vector normal to the fault plane (pointing upward) defined in the geographic reference system S
        this.normal_ = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)(); // constant values for each fault plane
        // (e_phi, e_theta) = unit vectors defining local reference frame tangent to the sphere in spherical coordinates
        this.e_phi = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.e_theta = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        // normalSp: unit vector normal to the fault plane (pointing upward) defined in the stress tensor reference system: S' = (X',Y',Z')=(s1,s3,s2)
        this.normalSp = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)(); // values should be recalculated for new stress tensors
        // (phiSp,thetaSp) : spherical coordinate angles defining the unit vector perpendicular to the fault plane (pointing upward in system S)
        //                 in the stress tensor reference system: S' = (X,Y,Z)
        // private phiSp:    number      // constant values for each fault plane
        // private thetaSp:  number      // values are recalculated for new stress tensors
        this.coordinatesSp = new _types_SphericalCoords__WEBPACK_IMPORTED_MODULE_1__.SphericalCoords();
        this.RTrot = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3)();
        // striation: unit vector pointing toward the measured striation in the geographic reference system: S = (X,Y,Z)
        // private striation:      Vector3 = newVector3D()      // constant value for each fault plane
        // striationSp: unit vector pointing toward the measured striation in the stress tensor reference system: S' = (X',Y',Z')
        this.striationSp = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)(); // values are recalculated for new stress tensors
        // stress: stress vector in the geographic reference system: S = (X,Y,Z)
        this.stress = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)(); // values are recalculated for new stress tensors
        // shearStressSp: shear stress vector in the geographic reference system: S = (X,Y,Z)
        this.shearStress = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.alphaStriaDeg = 0;
        this.alphaStria = 0;
        // e_striation_: unit vector in reference system S pointing toward the measured striation
        this.e_striation_ = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        // e_perp_striation_: unit vector in reference system S located in the fault plane and perpendicular to the measured striation
        this.e_perp_striation_ = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.newVector3D)();
        this.isUpLiftedBlock = false;
        this.strike = strike;
        this.dipDirection = dipDirection;
        this.dip = dip;
        this.faultSphericalCoords();
    }
    check({ displ, strain, stress }) {
        return stress !== undefined;
    }
    cost({ displ, strain, stress }) {
        return 1;
    }
    /**
     * General case to set the striation.
     *
     * Set the orientation of the striation in the fault plane, which can defined in two different ways (which are exclusive):
     * 1. Rake (or pitch) [0,90], measured from the strike direction, which points in one of the two opposite directions of the fault strike.
     *   Strike direction : (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
     * 2. For shallow-dipping planes (i.e., the compass inclinometer is inaccurate):
     *   Striae trend: [0, 360)
     */
    setStriation({ sensOfMovement, rake, strikeDirection }) {
        // check and set
        this.rake = rake;
        this.strikeDirection = strikeDirection;
        this.sensMouv = sensOfMovement;
        this.faultStriationAngle_A();
        this.faultStriationAngle_B();
        return this;
    }
    /**
     * Special case for shallow angle plane.
     */
    setStriationShallowPlane({ sensOfMovement, striationTrend }) {
        this.sensMouv = sensOfMovement;
        this.striationTrend = striationTrend;
        this.faultStriationAngle_A();
        this.faultStriationAngle_B();
        return this;
    }
    // TODO
    setStriationForVerticalPlane() {
    }
    // ------------------------------ PRIVATE
    faultSphericalCoords() {
        // Each fault is defined by a set of parameters as follows:
        //      The fault plane orientation is defined by three parameters:
        //      Fault strike: clockwise angle measured from the North direction [0, 360)
        //      Fault dip: [0, 90]
        //      Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
        // (phi,theta) : spherical coordinate angles defining the unit vector perpendicular to the fault plane (pointing upward)
        //                 in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
        // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
        // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
        //  Write functions relating trend and rake
        // The polar angle (or colatitude) theta is defined by the dip of the fault plane in radians:
        this.coordinates.theta = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.dip);
        // This function calculates the azimuth phi such that the right-handed local coordinate system in polar coordinates is located in the upper hemisphere.
        //      In other words, the radial unit vector is in the upper hemisphere.
        // The right-handed local reference system is specified by three unit vectors defined in the increasing radial, polar, and azimuthal directions (r, theta, and phi):
        //      The azimuthal angle phi is chosen in the direction of the fault dip (note that phi is different from the azimuth of the fault plane measured in the field) 
        //      The unit vector e_theta is parallel to the dip of the fault plane
        //      The unit vector e_phi is is parallel to the strike of the fault plane, and is oriented such that e_theta x e_phi = e_r (where x is the cross porduct )
        //      
        // The following 'if structure' calculates phi from the strike and dip direction of the fault plane:
        if (this.dip === 90) {
            // The fault plane is vertical and the dip direction is not defined
            if (this.strike <= 180) {
                // phi is in interval [0,PI]
                this.coordinates.phi = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else {
                // fault strike is in interval (PI,2 PI) and phi is in interval (PI,2 PI)
                this.coordinates.phi = 3 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
        }
        else if (this.strike === 0) { // The fault plane is not vertical and the dip direction is defined
            if (this.dipDirection === 0 /* Direction.E */) {
                this.coordinates.phi = 0;
            }
            else if (this.dipDirection === 1 /* Direction.W */) {
                this.coordinates.phi = Math.PI;
            }
            else {
                throw new Error(`dip direction is wrong. Should be E or W`);
            }
        }
        else if (this.strike < 90) {
            if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 5 /* Direction.SE */)) {
                // this.strike + this.coordinates.phi = 2Pi
                this.coordinates.phi = 2 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 7 /* Direction.NW */)) {
                // this.strike + this.coordinates.phi = Pi
                this.coordinates.phi = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else {
                throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
            }
        }
        else if (this.strike === 90) {
            if (this.dipDirection === 3 /* Direction.S */) {
                this.coordinates.phi = 3 * Math.PI / 2;
            }
            else if (this.dipDirection === 2 /* Direction.N */) {
                this.coordinates.phi = Math.PI / 2;
            }
            else {
                throw new Error(`dip direction is wrong. Should be N or S`);
            }
        }
        else if (this.strike < 180) {
            if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 6 /* Direction.SW */)) {
                // this.strike + this.coordinates.phi = 2Pi
                this.coordinates.phi = 2 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 4 /* Direction.NE */)) {
                // this.strike + this.coordinates.phi = Pi
                this.coordinates.phi = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else {
                throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
            }
        }
        else if (this.strike === 180) {
            if (this.dipDirection === 1 /* Direction.W */) {
                this.coordinates.phi = Math.PI;
            }
            else if (this.dipDirection === 0 /* Direction.E */) {
                this.coordinates.phi = 0;
            }
            else {
                throw new Error(`dip direction is wrong. Should be E or W`);
            }
        }
        else if (this.strike < 270) {
            if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 7 /* Direction.NW */)) {
                // this.strike + this.coordinates.phi = 2Pi
                this.coordinates.phi = 2 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 5 /* Direction.SE */)) {
                // this.strike + this.coordinates.phi = 3Pi
                this.coordinates.phi = 3 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else {
                throw new Error(`dip direction is wrong. Should be N, S, E, W, NW or SE`);
            }
        }
        else if (this.strike === 270) {
            if (this.dipDirection === 3 /* Direction.S */) {
                this.coordinates.phi = 3 * Math.PI / 2;
            }
            else if (this.dipDirection === 2 /* Direction.N */) {
                this.coordinates.phi = Math.PI / 2;
            }
            else {
                throw new Error(`dip direction is wrong. Should be N or S`);
            }
        }
        else if (this.strike < 360) {
            if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 4 /* Direction.NE */)) {
                // this.strike + this.coordinates.phi = 2Pi
                this.coordinates.phi = 2 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 6 /* Direction.SW */)) {
                // this.strike + this.coordinates.phi = 3Pi
                this.coordinates.phi = 3 * Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.strike);
            }
            else {
                throw new Error(`dip direction is wrong. Should be N, S, E, W, NE or SW`);
            }
        }
        else if (this.strike === 360) {
            if (this.dipDirection === 0 /* Direction.E */) {
                this.coordinates.phi = 0;
            }
            else if (this.dipDirection === 1 /* Direction.W */) {
                this.coordinates.phi = Math.PI;
            }
            else {
                throw new Error(`dip direction is wrong. Should be E or W`);
            }
        }
        else {
            throw new Error(`Strike is wrong. Should be in interval [0,360]`);
        }
        // The fault plane is defined by angles (phi, theta) in spherical coordinates.
        // normal: unit vector normal to the fault plane (pointing upward) defined in the geographic reference system: S = (X,Y,Z)
        this.normal_ = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.spherical2unitVectorCartesian)(this.coordinates);
        // --------------------------------------        
    }
    faultStriationAngle_A() {
        // Function calculating the striation angle in the local reference frame in polar coordinates from the rake
        //      The effect of fault movement on the striation is considered in function faultStriationAngle_B
        // Each fault is defined by a set of parameters as follows:
        //      The fault plane orientation is defined by three parameters:
        //      Fault strike: clockwise angle measured from the North direction [0, 360)
        //      Strike direction (optional): (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Fault dip: [0, 90]
        //      Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
        // The orientation of the striation in the fault plane can defined in two different ways (which are exclusive):
        // 1-   Rake (or pitch) [0,90], measured from the strike direction, which points in one of the two opposite directions of the fault strike.
        //      Strike direction : (N, E, S, W) or a combination of two direction (NE, SE, SW, NW).
        //      Note that the specified strike direction is used to determine the spatial orientation of the striation 
        // 2-   For shallow-dipping planes (i.e., the compass inclinometer is inaccurate):
        //      Striae trend: [0, 360)
        // alphaStria : striation angle measured in the local reference plane (e_phi, e_theta) indicating the motion of the outward block
        //      alphaStria is measured clockwise from e_phi, in interval [0, 2 PI) (this choice is consistent with the definition of the rake, which is measured from the fault strike)
        this.e_phi[0] = -Math.sin(this.coordinates.phi);
        this.e_phi[1] = Math.cos(this.coordinates.phi);
        this.e_phi[2] = 0;
        this.e_theta[0] = Math.cos(this.coordinates.theta) * Math.cos(this.coordinates.phi);
        this.e_theta[1] = Math.cos(this.coordinates.theta) * Math.sin(this.coordinates.phi);
        this.e_theta[2] = -Math.sin(this.coordinates.theta);
        // V[0] = Math.sin(this.coordinates.theta) * Math.cos( this.coordinates.phi )
        // V[1] = Math.sin(this.coordinates.theta) * Math.sin( this.coordinates.phi )
        // V[2] = Math.cos(this.coordinates.theta)
        // if structure for calculating the striation angle in the local reference frame in polar coordinates from the rake:
        // The following 'if structure' calculates phi from the strike and dip direction (if defined) of the fault plane:
        if (this.dip === 90) {
            // The fault plane is vertical and the dip direction is not defined
            if (this.strike === 0) {
                // phi = PI
                if (this.strikeDirection === 2 /* Direction.N */) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if (this.strikeDirection === 3 /* Direction.S */) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                }
            }
            else if (this.strike < 90) {
                // phi = PI - strike
                if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be  N, S, E, W, NE or SW`);
                }
            }
            else if (this.strike === 90) {
                // phi = PI/2
                if (this.strikeDirection === 0 /* Direction.E */) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if (this.strikeDirection === 1 /* Direction.W */) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                }
            }
            else if (this.strike < 180) {
                // phi = PI - strike
                if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SE or NW `);
                }
            }
            else if (this.strike === 180) {
                // phi = 0
                if (this.strikeDirection === 3 /* Direction.S */) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if (this.strikeDirection === 2 /* Direction.N */) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                }
            }
            else if (this.strike < 270) {
                // phi = 3 PI - strike
                if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SW or NE `);
                }
            }
            else if (this.strike === 270) {
                // phi = 3 PI / 2
                if (this.strikeDirection === 1 /* Direction.W */) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if (this.strikeDirection === 0 /* Direction.E */) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                }
            }
            else if (this.strike < 360) {
                // phi = 3 PI - strike
                if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, NW or SE `);
                }
            }
            else if (this.strike === 360) {
                // This case should not occur since in principle strike < 360
                // phi = PI
                if (this.strikeDirection === 2 /* Direction.N */) {
                    this.alphaStriaDeg = 180 - this.rake;
                    this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else if (this.strikeDirection === 3 /* Direction.S */) {
                    this.alphaStriaDeg = this.rake;
                    this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                }
                else {
                    throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                }
            }
            else {
                // This case should not occur since in principle strike < 360
                throw new Error(`fault strike is out of the expected interval [0,360)`);
            }
        }
        else { // The fault plane is not vertical and the dip direction is defined
            if (this.strike === 0) {
                if (this.dipDirection === 0 /* Direction.E */) {
                    if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = this.rake; // For testing the sense of mouvement of faults 
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = 180 - this.rake; // For testing the sense of mouvement of faults
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else if (this.dipDirection === 1 /* Direction.W */) {
                    if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be E or W`);
                }
            }
            else if (this.strike < 90) {
                if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 5 /* Direction.SE */)) {
                    if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, NE or SW `);
                    }
                }
                else if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 7 /* Direction.NW */)) {
                    if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be  N, S, E, W, NE or SW`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
                }
            }
            else if (this.strike === 90) {
                if (this.dipDirection === 3 /* Direction.S */) {
                    if (this.strikeDirection === 0 /* Direction.E */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 1 /* Direction.W */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                    }
                }
                else if (this.dipDirection === 2 /* Direction.N */) {
                    if (this.strikeDirection === 0 /* Direction.E */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 1 /* Direction.W */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N or S`);
                }
            }
            else if (this.strike < 180) {
                if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 6 /* Direction.SW */)) {
                    if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SE or NW `);
                    }
                }
                else if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 4 /* Direction.NE */)) {
                    if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SE or NW `);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N, S, E, W, SW or NE`);
                }
            }
            else if (this.strike === 180) {
                if (this.dipDirection === 1 /* Direction.W */) {
                    if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else if (this.dipDirection === 0 /* Direction.E */) {
                    if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be E or W`);
                }
            }
            else if (this.strike < 270) {
                if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 7 /* Direction.NW */)) {
                    if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SW or NE `);
                    }
                }
                else if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 5 /* Direction.SE */)) {
                    if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 6 /* Direction.SW */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 4 /* Direction.NE */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, SW or NE `);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
                }
            }
            else if (this.strike === 270) {
                if (this.dipDirection === 2 /* Direction.N */) {
                    if (this.strikeDirection === 1 /* Direction.W */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 0 /* Direction.E */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                    }
                }
                else if (this.dipDirection === 3 /* Direction.S */) {
                    if (this.strikeDirection === 1 /* Direction.W */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 0 /* Direction.E */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be E or W`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N or S`);
                }
            }
            else if (this.strike < 360) {
                if ((this.dipDirection === 2 /* Direction.N */) || (this.dipDirection === 0 /* Direction.E */) || (this.dipDirection === 4 /* Direction.NE */)) {
                    if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, NW or SE `);
                    }
                }
                else if ((this.dipDirection === 3 /* Direction.S */) || (this.dipDirection === 1 /* Direction.W */) || (this.dipDirection === 6 /* Direction.SW */)) {
                    if ((this.strikeDirection === 2 /* Direction.N */) || (this.strikeDirection === 1 /* Direction.W */) || (this.strikeDirection === 7 /* Direction.NW */)) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if ((this.strikeDirection === 3 /* Direction.S */) || (this.strikeDirection === 0 /* Direction.E */) || (this.strikeDirection === 5 /* Direction.SE */)) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N, S, E, W, NW or SE `);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be N, S, E, W, NE or SW`);
                }
            }
            else if (this.strike === 360) {
                // This case should not occur since in principle strike < 360
                if (this.dipDirection === 0 /* Direction.E */) {
                    if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else if (this.dipDirection === 1 /* Direction.W */) {
                    if (this.strikeDirection === 2 /* Direction.N */) {
                        this.alphaStriaDeg = 180 - this.rake;
                        this.alphaStria = Math.PI - (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else if (this.strikeDirection === 3 /* Direction.S */) {
                        this.alphaStriaDeg = this.rake;
                        this.alphaStria = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(this.rake);
                    }
                    else {
                        throw new Error(`Strike direction for measuring the rake is wrong. Should be N or S`);
                    }
                }
                else {
                    throw new Error(`dip direction is wrong. Should be E or W`);
                }
            }
            else {
                throw new Error(`fault strike is out of the expected interval [0,360)`);
            }
        }
    }
    faultStriationAngle_B() {
        // Function introducuing the effect of fault movement on the striation angle
        // This function is called after function faultStriationAngle_A
        // We calculate a unit vector e_striation pointing toward the measured striation
        // Sense of mouvement: For verification purposes, it is recommended to indicate both the dip-slip and strike-slip compoenents, when possible. 
        //      Dip-slip component:
        //          N = Normal fault, 
        //          I = Inverse fault or thrust
        //      Strike-slip componenet:
        //          RL = Right-Lateral fault
        //          LL = Left-Lateral fault
        // Sense of mouvement: N, I, RL, LL, N-RL, N-LL, I-RL, I-LL
        // this.alphaStriaDeg is in interval [0,180] according to function faultStriationAngle_A; 
        // This angle indicates the mouvement of the top (outward) block relative to the bottom (inner) block 
        // 'if structure' that modifies when required the striation angle according to the sense of mouvement of faults:
        if (this.dip === 90) {
            // The fault plane is vertical and only the strike-slip component of motion is defined
            // alphaStriaDeg calculated in function faultStriationAngle_B is in interval [0,PI]
            if ((this.alphaStriaDeg >= 0) && (this.alphaStriaDeg < 90)) {
                // alphaStriaDeg has a left-lateral strike-slip component 
                if (this.sensMouv === 3 /* SensOfMovement.RL */) {
                    // Fault movement is oriented opposite to the present value of the striation angle
                    this.alphaStriaDeg += 180;
                    this.alphaStria += Math.PI;
                }
                else if (this.sensMouv != 4 /* SensOfMovement.LL */) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be RL or LL`);
                }
            }
            else if (this.alphaStriaDeg === 90) { // Pure dip-slip mouvement
                // note that if alphaStriaDeg = 90 then the fault only has a dip-slip component and the direction of the uplifted block is requested 
                this.faultStriationUpliftedBlock();
            }
            else if (this.alphaStriaDeg <= 180) {
                // 90 < alphaStriaDeg <= 180 means that the fault is normal-right-lateral
                if (this.sensMouv === 4 /* SensOfMovement.LL */) {
                    // Fault movement is oriented opposite to the present value of the striation angle
                    this.alphaStriaDeg += 180;
                    this.alphaStria += Math.PI;
                }
                else if (this.sensMouv != 3 /* SensOfMovement.RL */) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be RL or LL`);
                }
            }
            else {
                throw new Error(`calculated striation alphaStriaDeg should be in interval [0,180]. Check routine faultStriationAngle_A`);
            }
        }
        else { // The fault plane is not vertical and both strike-slip and dip-slip components of motion are defined
            if (this.alphaStriaDeg === 0) { // Pure strike-slip mouvement
                // alphaStriaDeg = 0 means that the fault is left-lateral
                if (this.sensMouv === 3 /* SensOfMovement.RL */) {
                    // Fault movement is oriented opposite to the present value of the striation angle
                    this.alphaStriaDeg = 180; // Striation values are recalculated
                    this.alphaStria = Math.PI;
                }
                else if (this.sensMouv != 4 /* SensOfMovement.LL */) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be RL or LL`);
                }
            }
            else if (this.alphaStriaDeg < 90) { // Strike-slip and dip slip mouvement
                // 0 < alphaStriaDeg < 90 means that the fault is normal-left-lateral
                if ((this.sensMouv === 3 /* SensOfMovement.RL */) || (this.sensMouv === 2 /* SensOfMovement.I */) || (this.sensMouv === 7 /* SensOfMovement.I_RL */)) {
                    this.alphaStriaDeg += 180;
                    this.alphaStria += Math.PI;
                }
                else if ((this.sensMouv !== 4 /* SensOfMovement.LL */) && (this.sensMouv !== 1 /* SensOfMovement.N */) && (this.sensMouv !== 6 /* SensOfMovement.N_LL */)) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be LL or N or N-LL or RL or I or I-RL`);
                }
            }
            else if (this.alphaStriaDeg === 90) { // Pure dip-slip mouvement
                // alphaStriaDeg = 90 means that the fault is normal
                if (this.sensMouv === 2 /* SensOfMovement.I */) {
                    // Fault movement is oriented opposite to the present value of the striation angle
                    this.alphaStriaDeg = 270; // Striation values are recalculated
                    this.alphaStria = 3 * Math.PI / 2;
                }
                else if (this.sensMouv != 1 /* SensOfMovement.N */) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be N or I`);
                }
            }
            else if (this.alphaStriaDeg < 180) { // Strike-slip and dip slip mouvement
                // 90 < alphaStriaDeg < 180 means that the fault is normal-right-lateral
                if ((this.sensMouv === 4 /* SensOfMovement.LL */) || (this.sensMouv === 2 /* SensOfMovement.I */) || (this.sensMouv === 8 /* SensOfMovement.I_LL */)) {
                    this.alphaStriaDeg += 180;
                    this.alphaStria += Math.PI;
                }
                else if ((this.sensMouv != 3 /* SensOfMovement.RL */) && (this.sensMouv != 1 /* SensOfMovement.N */) && (this.sensMouv === 5 /* SensOfMovement.N_RL */)) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be LL or I or I-LL or RL or N or N-RL`);
                }
            }
            else if (this.alphaStriaDeg === 180) { // Pure strike-slip mouvement
                // alphaStriaDeg = 180 means that the fault is right-lateral
                if (this.sensMouv === 4 /* SensOfMovement.LL */) {
                    // Fault movement is oriented opposite to the present value of the striation angle
                    this.alphaStriaDeg = 0; // Striation values are recalculated
                    this.alphaStria = 0;
                }
                else if (this.sensMouv != 3 /* SensOfMovement.RL */) {
                    throw new Error(`sense of mouvement is not consistent with fault data. Should be RL or LL`);
                }
            }
            else {
                throw new Error(`calculated striation alphaStriaDeg should be in interval [0,180]. Check routine faultStriationAngle_A`);
            }
        }
        // Calculate in reference system S the unit vector e_striation pointing toward the measured striation BB
        this.e_striation_[0] = Math.cos(this.alphaStria) * this.e_phi[0] + Math.sin(this.alphaStria) * this.e_theta[0];
        this.e_striation_[1] = Math.cos(this.alphaStria) * this.e_phi[1] + Math.sin(this.alphaStria) * this.e_theta[1];
        this.e_striation_[2] = Math.cos(this.alphaStria) * this.e_phi[2] + Math.sin(this.alphaStria) * this.e_theta[2];
        // Calculate in reference system S the unit vector e_perp_striation_ located on the fault plane and perpendiculat to the striation.
        // This vector is necessary for calculating the misfit angle for criterai involving friction.
        // The local coord system (e_striation_, e_perp_striation_, normal) is right handed
        this.e_perp_striation_ = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.crossProduct)({ U: this.normal, V: this.e_striation_ });
    }
    faultStriationUpliftedBlock() {
        // The fault plane is vertical and the rake is 90°: this.dip = 90, this.rake = 90, this.alphaStriaDeg = 90
        //Thus, the striation is defined by an additional parameter: 
        // To calculate the orientation of the striation the user must indicate for example the direction of the uplifted block:
        //      upLiftedBlock: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW)
        if (this.strike === 0) {
            if (this.upliftedBlock === 1 /* Direction.W */) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if (this.upliftedBlock !== 0 /* Direction.E */) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be E or W`);
            }
        }
        else if (this.strike < 90) {
            if ((this.upliftedBlock === 2 /* Direction.N */) || (this.upliftedBlock === 1 /* Direction.W */) || (this.upliftedBlock === 7 /* Direction.NW */)) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if ((this.upliftedBlock !== 3 /* Direction.S */) && (this.upliftedBlock !== 0 /* Direction.E */) && (this.upliftedBlock !== 5 /* Direction.SE */)) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N, S, E, W, SE or NW`);
            }
        }
        else if (this.strike === 90) {
            if (this.upliftedBlock === 2 /* Direction.N */) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if (this.upliftedBlock !== 3 /* Direction.S */) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N or S`);
            }
        }
        else if (this.strike < 180) {
            if ((this.upliftedBlock === 2 /* Direction.N */) || (this.upliftedBlock === 0 /* Direction.E */) || (this.upliftedBlock === 4 /* Direction.NE */)) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if ((this.upliftedBlock !== 3 /* Direction.S */) && (this.upliftedBlock !== 1 /* Direction.W */) && (this.upliftedBlock !== 6 /* Direction.SW */)) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N, S, E, W, NE or SW`);
            }
        }
        else if (this.strike === 180) {
            if (this.upliftedBlock === 0 /* Direction.E */) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if (this.upliftedBlock !== 1 /* Direction.W */) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be E or W`);
            }
        }
        else if (this.strike < 270) {
            if ((this.upliftedBlock === 3 /* Direction.S */) || (this.upliftedBlock === 0 /* Direction.E */) || (this.upliftedBlock === 5 /* Direction.SE */)) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if ((this.upliftedBlock !== 2 /* Direction.N */) && (this.upliftedBlock !== 1 /* Direction.W */) && (this.upliftedBlock !== 7 /* Direction.NW */)) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N, S, E, W, SE or NW`);
            }
        }
        else if (this.strike === 270) {
            if (this.upliftedBlock === 3 /* Direction.S */) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if (this.upliftedBlock !== 2 /* Direction.N */) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N or S`);
            }
        }
        else if (this.strike < 360) {
            if ((this.upliftedBlock === 3 /* Direction.S */) || (this.upliftedBlock === 1 /* Direction.W */) || (this.upliftedBlock === 6 /* Direction.SW */)) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if ((this.upliftedBlock !== 2 /* Direction.N */) && (this.upliftedBlock !== 0 /* Direction.E */) && (this.upliftedBlock !== 4 /* Direction.NE */)) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be N, S, E, W, NE or SW`);
            }
        }
        else if (this.strike === 360) {
            if (this.upliftedBlock === 1 /* Direction.W */) {
                this.alphaStriaDeg = 270;
                this.alphaStria = 3 * Math.PI / 2;
            }
            else if (this.upliftedBlock !== 0 /* Direction.E */) {
                throw new Error(`The orientation of the uplifted block is wrong. Should be E or W`);
            }
        }
        else {
            throw new Error(`fault strike is out of the expected interval [0,360)`);
        }
    }
    createUpLiftedBlock() {
        const f = new Fault({ strike: 0, dipDirection: 0 /* Direction.E */, dip: 90 }); // TODO: params in ctor
        f.setStriation({ strikeDirection: 2 /* Direction.N */, rake: 90, sensOfMovement: 9 /* SensOfMovement.UKN */ });
        /**
        * There is one particular case in which the sens of mouvement has to be defined with a different parameter:
        * A vertical plane with rake = 90.
        * In such situation the user must indicate for example the direction of the uplifted block:
        *      upLiftedBlock: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
        */
        return f;
    }
    faultSphericalCoordsSP() {
        // Calculate the spherical coordinates of the unit vector normal to the plane in reference system S'
        // normalSp: unit vector normal to the fault plane (pointing upward) defined in the stress tensor reference system: S' = (X',Y',Z')=(sigma 1,sigma 3,sigma 2)
        //          values should be recalculated for new stress tensors    
        // Let Rrot be the rotation tensor R between reference systems S and S', such that:
        //      V' = R V,  where V and V' are the same vector defined in reference frames S and S', respectively
        this.normalSp = (0,_types_math__WEBPACK_IMPORTED_MODULE_0__.tensor_x_Vector)({ T: this.RTrot, V: this.normal });
        if (this.normalSp[0] > 0) {
            if (this.normalSp[1] >= 0) {
                // phiSp is in interval [0, Pi/2)
                this.coordinatesSp.phi = Math.atan(this.normalSp[1] / this.normalSp[0]);
            }
            else {
                // phiSp is in interval (3Pi/2, 2Pi)
                // atan is probably defined in interval (-Pi/2, Pi/2)
                this.coordinatesSp.phi = 2 * Math.PI + Math.atan(this.normalSp[1] / this.normalSp[0]);
            }
        }
        else if (this.normalSp[0] < 0) {
            if (this.normalSp[1] >= 0) {
                // phiSp is in interval (Pi/2, Pi]
                this.coordinatesSp.phi = Math.atan(this.normalSp[1] / this.normalSp[0]) + Math.PI;
            }
            else {
                // phiSp is defined in interval [Pi, 3Pi/2)
                this.coordinatesSp.phi = Math.atan(this.normalSp[1] / this.normalSp[0]) + Math.PI;
            }
        }
        else {
            if (this.normalSp[1] > 0) {
                // phiSp = Pi/2
                this.coordinatesSp.phi = Math.PI / 2;
            }
            else {
                // phiSp = 3Pi/2
                this.coordinatesSp.phi = 3 * Math.PI / 2;
            }
        }
    }
    faultNormalVectorSp() {
        /**
         *  (phiSp,thetaSp) : spherical coordinate angles defining the unit vector perpendicular to the fault plane (pointing upward in system S)
         *               in the stress tensor reference system: S' = (X,Y,Z)
         *  These angles are recalculated from the new stress tensors
         */
    }
    /**
     * Rotate the tensor about an angle...
     * @param rotAx_phi
     */
    vector_rotation(rotAx_phi) {
        // this.x = Math.sin( rotAx_phi);
    }
}
// Stress parameters:
// For a given stress tensor, we calculate the stress components:
// Step 1:
//      n' = R n,  where n and n' are vectors in reference frames S and S'
/*
this.normalSp[0] = R[0,0] * this.normal[0] + R[0,1] * this.normal[1] + R[0,2] * this.normal[2]
this.normalSp[1] = R[1,0] * this.normal[0] + R[1,1] * this.normal[1] + R[1,2] * this.normal[2]
this.normalSp[2] = R[2,0] * this.normal[0] + R[2,1] * this.normal[1] + R[2,2] * this.normal[2]

// Step 2:
// The principal stress values are defined according to the rock mechanics sign convention (positive values for compressive stresses)
const sigma_1 = - this.lambda[0]    // Principal stress in X direction
const sigma_2 = - this.lambda[2]    // Principal stress in Z direction
const sigma_3 = - this.lambda[1]    // Principal stress in Y direction

// Calculate the normal and shear stress components of the fault plane using coordinates in reference system S':
// S' = (X',Y',Z') is the principal stress reference frame, parallel to (sigma_1, sigma_3, sigma_2);
// The stress magnitude is obtained from the sum of the squared components
let this.Stress = Math.sqrt( sigma_1**2 * np[0]**2 + sigma_3**2 * np[1]**2 + sigma_2**2 * np[2]**2 )
// The signed normal stress is obtatined form the scalar product of the normal and stress vectors
// The normal stress is positive since (s1,s2,s3) = (1,R,0)
let this.normalStress = sigma_1 * np[0]**2 + sigma_3 * np[1]**2 + sigma_2 * np[2]**2
// The shear stress
let this.shearStress = sigma_1 * np[0]**2 + sigma_3 * np[1]**2 + sigma_2 * np[2]**2
*/


/***/ }),

/***/ "./lib/utils/fromAnglesToNormal.ts":
/*!*****************************************!*\
  !*** ./lib/utils/fromAnglesToNormal.ts ***!
  \*****************************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "fromAnglesToNormal": () => (/* binding */ fromAnglesToNormal)
/* harmony export */ });
/* harmony import */ var _types__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types */ "./lib/types/index.ts");

function fromAnglesToNormal({ strike, dip, dipDirection }) {
    // Each fault is defined by a set of parameters as follows:
    //      The fault plane orientation is defined by three parameters:
    //      Fault strike: clockwise angle measured from the North direction [0, 360)
    //      Fault dip: [0, 90]
    //      Dip direction: (N, E, S, W) or a combination of two directions (NE, SE, SW, NW).
    // (phi,theta) : spherical coordinate angles defining the unit vector perpendicular to the fault plane (pointing upward)
    //                 in the geographic reference system: S = (X,Y,Z) = (E,N,Up)
    // phi : azimuthal angle in interval [0, 2 PI), measured anticlockwise from the X axis (East direction) in reference system S
    // theta: colatitude or polar angle in interval [0, PI/2], measured downward from the zenith (upward direction)
    //  Write functions relating trend and rake
    // The polar angle (or colatitude) theta is defined by the dip of the fault plane in radians:
    const coordinates = new _types__WEBPACK_IMPORTED_MODULE_0__.SphericalCoords;
    coordinates.theta = (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(dip);
    // This function calculates the azimuth phi such that the right-handed local coordinate system in polar coordinates is located in the upper hemisphere.
    //      In other words, the radial unit vector is in the upper hemisphere.
    // The right-handed local reference system is specified by three unit vectors defined in the increasing radial, polar, and azimuthal directions (r, theta, and phi):
    //      The azimuthal angle phi is chosen in the direction of the fault dip (note that phi is different from the azimuth of the fault plane measured in the field) 
    //      The unit vector e_theta is parallel to the dip of the fault plane
    //      The unit vector e_phi is is parallel to the strike of the fault plane, and is oriented such that e_theta x e_phi = e_r (where x is the cross porduct )
    //      
    // The following 'if structure' calculates phi from the strike and dip direction of the fault plane:
    if (dip === 90) {
        // The fault plane is vertical and the dip direction is not defined
        if (strike <= 180) {
            // phi is in interval [0,PI]
            coordinates.phi = Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else {
            // fault strike is in interval (PI,2 PI) and phi is in interval (PI,2 PI)
            coordinates.phi = 3 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
    }
    else if (strike === 0) { // The fault plane is not vertical and the dip direction is defined
        if (dipDirection === 0 /* Direction.E */) {
            coordinates.phi = 0;
        }
        else if (dipDirection === 1 /* Direction.W */) {
            coordinates.phi = Math.PI;
        }
        else {
            throw new Error(`dip direction is wrong. Should be E or W`);
        }
    }
    else if (strike < 90) {
        if ((dipDirection === 3 /* Direction.S */) || (dipDirection === 0 /* Direction.E */) || (dipDirection === 5 /* Direction.SE */)) {
            // strike + coordinates.phi = 2Pi
            coordinates.phi = 2 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else if ((dipDirection === 2 /* Direction.N */) || (dipDirection === 1 /* Direction.W */) || (dipDirection === 7 /* Direction.NW */)) {
            // strike + coordinates.phi = Pi
            coordinates.phi = Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else {
            throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
        }
    }
    else if (strike === 90) {
        if (dipDirection === 3 /* Direction.S */) {
            coordinates.phi = 3 * Math.PI / 2;
        }
        else if (dipDirection === 2 /* Direction.N */) {
            coordinates.phi = Math.PI / 2;
        }
        else {
            throw new Error(`dip direction is wrong. Should be N or S`);
        }
    }
    else if (strike < 180) {
        if ((dipDirection === 3 /* Direction.S */) || (dipDirection === 1 /* Direction.W */) || (dipDirection === 6 /* Direction.SW */)) {
            // strike + coordinates.phi = 2Pi
            coordinates.phi = 2 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else if ((dipDirection === 2 /* Direction.N */) || (dipDirection === 0 /* Direction.E */) || (dipDirection === 4 /* Direction.NE */)) {
            // strike + coordinates.phi = Pi
            coordinates.phi = Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else {
            throw new Error(`dip direction is wrong. Should be N, S, E, W, SE or NW`);
        }
    }
    else if (strike === 180) {
        if (dipDirection === 1 /* Direction.W */) {
            coordinates.phi = Math.PI;
        }
        else if (dipDirection === 0 /* Direction.E */) {
            coordinates.phi = 0;
        }
        else {
            throw new Error(`dip direction is wrong. Should be E or W`);
        }
    }
    else if (strike < 270) {
        if ((dipDirection === 2 /* Direction.N */) || (dipDirection === 1 /* Direction.W */) || (dipDirection === 7 /* Direction.NW */)) {
            // strike + coordinates.phi = 2Pi
            coordinates.phi = 2 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else if ((dipDirection === 3 /* Direction.S */) || (dipDirection === 0 /* Direction.E */) || (dipDirection === 5 /* Direction.SE */)) {
            // strike + coordinates.phi = 3Pi
            coordinates.phi = 3 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else {
            throw new Error(`dip direction is wrong. Should be N, S, E, W, NW or SE`);
        }
    }
    else if (strike === 270) {
        if (dipDirection === 3 /* Direction.S */) {
            coordinates.phi = 3 * Math.PI / 2;
        }
        else if (dipDirection === 2 /* Direction.N */) {
            coordinates.phi = Math.PI / 2;
        }
        else {
            throw new Error(`dip direction is wrong. Should be N or S`);
        }
    }
    else if (strike < 360) {
        if ((dipDirection === 2 /* Direction.N */) || (dipDirection === 0 /* Direction.E */) || (dipDirection === 4 /* Direction.NE */)) {
            // strike + coordinates.phi = 2Pi
            coordinates.phi = 2 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else if ((dipDirection === 3 /* Direction.S */) || (dipDirection === 1 /* Direction.W */) || (dipDirection === 6 /* Direction.SW */)) {
            // strike + coordinates.phi = 3Pi
            coordinates.phi = 3 * Math.PI - (0,_types__WEBPACK_IMPORTED_MODULE_0__.deg2rad)(strike);
        }
        else {
            throw new Error(`dip direction is wrong. Should be N, S, E, W, NE or SW`);
        }
    }
    else if (strike === 360) {
        if (dipDirection === 0 /* Direction.E */) {
            coordinates.phi = 0;
        }
        else if (dipDirection === 1 /* Direction.W */) {
            coordinates.phi = Math.PI;
        }
        else {
            throw new Error(`dip direction is wrong. Should be E or W`);
        }
    }
    else {
        throw new Error(`Strike is wrong. Should be in interval [0,360]`);
    }
    // The fault plane is defined by angles (phi, theta) in spherical coordinates.
    // normal: unit vector normal to the fault plane (pointing upward) defined in the geographic reference system: S = (X,Y,Z)
    return (0,_types__WEBPACK_IMPORTED_MODULE_0__.spherical2unitVectorCartesian)(coordinates);
}


/***/ }),

/***/ "./lib/utils/index.ts":
/*!****************************!*\
  !*** ./lib/utils/index.ts ***!
  \****************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "Direction": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.Direction),
/* harmony export */   "Fault": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.Fault),
/* harmony export */   "SensOfMovement": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.SensOfMovement),
/* harmony export */   "directionExists": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.directionExists),
/* harmony export */   "faultParams": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.faultParams),
/* harmony export */   "getDirectionFromString": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.getDirectionFromString),
/* harmony export */   "getSensOfMovementFromString": () => (/* reexport safe */ _Fault__WEBPACK_IMPORTED_MODULE_0__.getSensOfMovementFromString),
/* harmony export */   "trimAll": () => (/* reexport safe */ _trimAlls__WEBPACK_IMPORTED_MODULE_1__.trimAll)
/* harmony export */ });
/* harmony import */ var _Fault__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Fault */ "./lib/utils/Fault.ts");
/* harmony import */ var _trimAlls__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./trimAlls */ "./lib/utils/trimAlls.ts");




/***/ }),

/***/ "./lib/utils/trimAlls.ts":
/*!*******************************!*\
  !*** ./lib/utils/trimAlls.ts ***!
  \*******************************/
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "trimAll": () => (/* binding */ trimAll)
/* harmony export */ });
function trimAll(s) {
    return s
        .replace(/\s+/g, ' ')
        .replace(/^\s+|\s+$/, '')
        .replace('\t', ' ')
        .trimEnd();
}


/***/ }),

/***/ "@youwol/math":
/*!*******************************!*\
  !*** external "@youwol/math" ***!
  \*******************************/
/***/ ((module) => {

module.exports = __WEBPACK_EXTERNAL_MODULE__youwol_math__;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/compat get default export */
/******/ 	(() => {
/******/ 		// getDefaultExport function for compatibility with non-harmony modules
/******/ 		__webpack_require__.n = (module) => {
/******/ 			var getter = module && module.__esModule ?
/******/ 				() => (module['default']) :
/******/ 				() => (module);
/******/ 			__webpack_require__.d(getter, { a: getter });
/******/ 			return getter;
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/define property getters */
/******/ 	(() => {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = (exports, definition) => {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	(() => {
/******/ 		__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
/******/ 	})();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	(() => {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = (exports) => {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	})();
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
// This entry need to be wrapped in an IIFE because it need to be isolated against other modules in the chunk.
(() => {
/*!******************!*\
  !*** ./index.ts ***!
  \******************/
__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "CompactionShearBandsKin": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.CompactionShearBandsKin),
/* harmony export */   "ConjugateFaultsKin": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.ConjugateFaultsKin),
/* harmony export */   "Curve3D": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.Curve3D),
/* harmony export */   "Data": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.Data),
/* harmony export */   "DataDescription": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.DataDescription),
/* harmony export */   "DataFactory": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.DataFactory),
/* harmony export */   "Direction": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.Direction),
/* harmony export */   "EquipotentialCurve": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.EquipotentialCurve),
/* harmony export */   "ExtensionFracture": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.ExtensionFracture),
/* harmony export */   "Fault": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.Fault),
/* harmony export */   "FractureStrategy": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.FractureStrategy),
/* harmony export */   "GridSearch": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.GridSearch),
/* harmony export */   "IntegralCurve": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.IntegralCurve),
/* harmony export */   "InverseMethod": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.InverseMethod),
/* harmony export */   "MasterStress": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.MasterStress),
/* harmony export */   "MohrCoulombCurve": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.MohrCoulombCurve),
/* harmony export */   "MonteCarlo": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.MonteCarlo),
/* harmony export */   "PoleCoords": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.PoleCoords),
/* harmony export */   "SearchMethodFactory": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.SearchMethodFactory),
/* harmony export */   "SensOfMovement": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.SensOfMovement),
/* harmony export */   "SphericalCoords": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.SphericalCoords),
/* harmony export */   "StressTensor": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.StressTensor),
/* harmony export */   "StriatedPlaneFriction1": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.StriatedPlaneFriction1),
/* harmony export */   "StriatedPlaneFriction2": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.StriatedPlaneFriction2),
/* harmony export */   "StriatedPlaneKin": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.StriatedPlaneKin),
/* harmony export */   "StriatedPlaneProblemType": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.StriatedPlaneProblemType),
/* harmony export */   "add_Vectors": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.add_Vectors),
/* harmony export */   "angularDifStriations": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.angularDifStriations),
/* harmony export */   "arcCircle": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.arcCircle),
/* harmony export */   "cloneMatrix3x3": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.cloneMatrix3x3),
/* harmony export */   "cloneMisfitCriteriunSolution": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.cloneMisfitCriteriunSolution),
/* harmony export */   "constant_x_Vector": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.constant_x_Vector),
/* harmony export */   "createDefaultSolution": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.createDefaultSolution),
/* harmony export */   "crossProduct": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.crossProduct),
/* harmony export */   "decodeCSV": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.decodeCSV),
/* harmony export */   "deg2rad": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.deg2rad),
/* harmony export */   "directionExists": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.directionExists),
/* harmony export */   "faultParams": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.faultParams),
/* harmony export */   "faultStressComponents": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.faultStressComponents),
/* harmony export */   "getDirectionFromString": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.getDirectionFromString),
/* harmony export */   "getSensOfMovementFromString": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.getSensOfMovementFromString),
/* harmony export */   "lineSphericalCoords": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.lineSphericalCoords),
/* harmony export */   "minRotAngleRotationTensor": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.minRotAngleRotationTensor),
/* harmony export */   "mohrCircleLine": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.mohrCircleLine),
/* harmony export */   "mohrCirclePoint": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.mohrCirclePoint),
/* harmony export */   "multiplyTensors": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.multiplyTensors),
/* harmony export */   "newMatrix3x3": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3),
/* harmony export */   "newMatrix3x3Identity": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.newMatrix3x3Identity),
/* harmony export */   "newPoint3D": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.newPoint3D),
/* harmony export */   "newVector3D": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.newVector3D),
/* harmony export */   "normalVector": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.normalVector),
/* harmony export */   "normalizeVector": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.normalizeVector),
/* harmony export */   "normalizedCrossProduct": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.normalizedCrossProduct),
/* harmony export */   "properRotationTensor": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.properRotationTensor),
/* harmony export */   "rad2deg": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.rad2deg),
/* harmony export */   "rotationParamsFromRotTensor": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.rotationParamsFromRotTensor),
/* harmony export */   "scalarProduct": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.scalarProduct),
/* harmony export */   "scalarProductUnitVectors": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.scalarProductUnitVectors),
/* harmony export */   "setValueInUnitInterval": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.setValueInUnitInterval),
/* harmony export */   "signedAngularDifStriations": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.signedAngularDifStriations),
/* harmony export */   "spherical2unitVectorCartesian": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.spherical2unitVectorCartesian),
/* harmony export */   "stressTensorDelta": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.stressTensorDelta),
/* harmony export */   "stressTensorPrincipalAxes": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.stressTensorPrincipalAxes),
/* harmony export */   "tensor_x_Vector": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.tensor_x_Vector),
/* harmony export */   "transposeTensor": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.transposeTensor),
/* harmony export */   "trend2phi": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.trend2phi),
/* harmony export */   "trimAll": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.trimAll),
/* harmony export */   "unitVectorCartesian2Spherical": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.unitVectorCartesian2Spherical),
/* harmony export */   "vectorMagnitude": () => (/* reexport safe */ _lib__WEBPACK_IMPORTED_MODULE_0__.vectorMagnitude)
/* harmony export */ });
/* harmony import */ var _lib__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./lib */ "./lib/index.ts");


})();

/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=stress.js.map