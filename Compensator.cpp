#include "Compensator.h"

Compensator::Compensator(unsigned char _mode, double _r, double _s, double _a, double _t, double _method = Forward)
{
	this->mode = _mode;
	this->method = _method;

	this->a = _a;
	this->r = _r;
	this->s = _s;
	this->t = _t;

	if (_mode != LeadLag) this->s = this->r;
	this->st = this->s * this->t;
	this->rt = this->r * this->t;
	this->art = this->a * this->rt;
	this->stpa = this->st / this->a;
}

double Compensator::update(double ek)
{
	double ck = 0;
	switch (this->mode)
	{
	case Lead:
		switch (this->method)
		{
		case Forward:
			// ck = (1+rt)*ckm1 + ek - (1+art)*ekm1
			ck = (1 + this->rt) * this->ckm1 + ek - (1 + this->art) * this->ekm1;
			break;
		case Backward:
			// ck = (ckm1 + (1-art)*ek -ekm1)/(1-rt)
			ck = (this->ckm1 + (1 - this->art) * ek - this->ekm1) / (1 - this->rt);
			break;
		case Twisten:
			// ck = ((1+art/2)*ckm1 + (1-art/2)*ek - (1+art/2)*ekm1 )/(1-rt/2)
			ck = ((1 + this->art / 2) * this->ckm1 + (1 - 0.5 * this->art) * ek - (1 + 0.5 * this->art)
				* this->ekm1) / (1 - 0.5 * this->rt);
			break;
		default:
			// Forward euler again
			// ck = (1+rt)*ckm1 + ek - (1+art)*ekm1
			ck = (1 + this->rt) * this->ckm1 + ek - (1 + this->art) * this->ekm1;
			break;
		};
		break;

	case Lag:
		switch (this->method)
		{
		case Forward:
			// ck = (1-rt)ckm1 + (ek - (1-art)*ekm1)/a
			ck = (1 - this->rt) * this->ckm1 + (ek - (1 - this->art) * this->ekm1) / this->a;
			break;
		case Backward:
			// ck = (ckm1 + ( (1+art/2)*ek -ekm1 )/a )/(1+rt)
			ck = (this->ckm1 + ((1 + this->art / 2) * ek - this->ekm1) / this->a) / (1 + this->rt);
			break;
		case Twisten:
			// ck = ( (1-rt/2)*ckm1 + ( (1+art/2)*ek - (1-art/2)*ekm1 )/a )/(1+rt/2)
			ck = ((1 - 0.5 * this->rt) * this->ckm1 + ((1 + 0.5 * this->art) * ek - (1 - 0.5 * this->art)
				* this->ekm1) / this->a) / (1 + 0.5 * this->rt);
			break;
		default:
			// Forward euler again
			// ck = (1-rt)ckm1 + (ek - (1-art)*ekm1)/a
			ck = (1 - this->rt) * this->ckm1 + (ek - (1 - this->art) * this->ekm1) / this->a;
			break;
		};
		break;

	case LeadLag:


		switch (this->method)
		{
		case Forward:
			// ck = 
			ck = (this->st - 1) * (this->rt - 1) * this->ckm2 + (this->rt + this->st - 2) * this->ckm1;
			ck = ek + (this->art - 1) * (this->stpa - 1) * this->ekm2 + (this->art + this->stpa - 2) * this->ekm1 - ck;
			break;
		case Backward:
			// ck = 
			ck = (this->st + this->rt + 2) * this->ckm1 - this->ckm2 + (this->art + 1) * (this->st / a + 1) * ek
				- (this->art + this->st / a + 2) * this->ekm1 + this->ekm2;
			ck = ck / ((this->st + 1) * (this->rt + 1));
			break;
		case Twisten:
			// ck = 
			ck = ((0.5 * this->rt - 1) * (0.5 * this->st + 1) + (0.5 * this->st - 1) * (0.5 * this->rt + 1)) * this->ckm1;
			ck = ck + (0.5 * this->rt - 1) * (0.5 * this->st - 1) * this->ckm2;
			ck = -ck + (0.25 * this->st * this->rt + 0.5 * (this->art + this->stpa) + 1) * ek + (0.5 * this->st * this->rt - 2) * this->ekm1;
			ck = ck + (0.5 * this->art - 1) * (0.5 * this->stpa - 1) * this->ekm2;
			ck = ck / ((0.5 * this->rt + 1) * (0.5 * this->st + 1));
			break;
		default:
			// Forward euler again
			// ck = 
			ck = (this->st - 1) * (this->rt - 1) * this->ckm2 + (this->rt + this->st - 2) * this->ckm1;
			ck = ek + (this->art - 1) * (this->stpa - 1) * this->ekm2 + (this->art + this->stpa - 2) * this->ekm1 - ck;
			break;
		};
		break;
	};

	// update ckm1 and ekm1 for next calculation
	this->ckm2 = this->ckm1;
	this->ekm2 = this->ekm1;
	this->ckm1 = ck;
	this->ekm1 = ek;

	return ck;
}