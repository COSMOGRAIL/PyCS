Warning & Disclaimer
====================


By *itself*, the use of PyCS does not guarantee realistic time-delay uncertainty estimates. All we have (blindly) demonstrated regarding the quality of PyCS uncertainty estimates assumes a careful generation of the simulated light curves. Keep in mind that the PyCS approach involves setting values or ranges for several parameters. In principle, the impact of these settings has to be assessed on each light curve set, as described in the PyCS paper (`Tewes et al. 2013 <http://dx.doi.org/10.1051/0004-6361/201220123>`_).

.. warning:: In particular, obtaining uncertainty estimates by running on simulated light curves that all have a same true delay (i.e., setting ``truetsr = 0``, in PyCS terminology) can lead to **vastly underestimated error bars**.

It is of **high importance** to use simulations with a range of plausible true time delays to tune and/or verify the accuracy and precision of time-delay estimators. Tests on simulations with only a single true time delay do not probe reliably the quality of a time-delay estimation, as many time-delay estimators are prone to responding unsteadily to the true delay. Do not hesitate to contact us (see front page of this documentation) in case of doubts!

.. note:: Please try to avoid publishing time delays with overly optimistic uncertainty estimates. Seemingly accurate measurements might be propagated into unrealistic cosmological inferences, and harm the community. Be critical, and aim for the most comprehensive instead of the smallest error bar.

Thanks!