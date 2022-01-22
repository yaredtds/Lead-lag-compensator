#pragma once

#define Lead 0
#define Lag 1
#define LeadLag 2

#define Forward 0
#define Backward 1
#define Twisten 2

class Compensator
{
private:
	double a; // coefficient of compensator
	double r; // coefficient of compensator
	double s; // coefficient of lag compensator @ lead-lag
	double t; // sampling time
	unsigned char method; // Integral approx
	unsigned char mode;	// compensator mode of operation

	double ckm1 = 0; // control effort   @t-1
	double ekm1 = 0; // exitation signal @t-1
	double ekm2 = 0;
	double ckm2 = 0;

	double stpa = 0;
	double art = 0;  // product of a r and t
	double rt = 0;   // product of rt
	double st = 0;
public:
	void set_a(double _a = 0) { this->a = _a; }
	void set_r(double _r = 0) { this->r = _r; }
	void set_s(double _s = 0) { this->r = _s; }
	void set_t(double _t = 1) { this->t = _t; }
	void set_mode(char _mode) { this->mode = _mode; }
	void set_method(char _method) { this->method = _method; }

	double get_a() { return this->a; }
	double get_r() { return this->r; }
	double get_s() { return this->s; }
	double get_t() { return this->t; }

	Compensator(unsigned char _mode, double _r, double _s, double _a, double _t, double _method);
	double update(double ek);
};

