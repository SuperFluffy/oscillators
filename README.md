Coupled oscillators
===================

Networks of coupled oscillators are an object of active research in a massive
number of different scientific fields, among them statistical physics, biology,
cognitive science, and mathematics.

At the same time there is an odd absence of educational codes out there,
showing how to simulate such systems — shoot me a line if you think I have
missed a good source.

While learning about coupled oscillator networks, I will try and put together
some simple codes to simulate them. This repository is meant to organize my own
workflow, and I hope that somebody looking into coupled oscillator networks
might find them useful as well.

Types of oscillators:
* Kuramoto
* Wilson-Cowan
* van-der-Poel
* Spiking
  * Integrate-and-fire
  * Hodgkin-Huxley
  * FitzHugh-Nagumo


Kuramo model
------------

The Kuramoto model is a model of coupled phase oscillators. It is a hugely
popular model, while being quite easy to understand for a newbie in the field
(take this with a grain of salt: the equation itself looks quite straight
forward).

The implementation is based on [this tutorial](http://tutorials.siam.org/dsweb/cotutorial/).

