function [neuronfire,V,P]=neuron(w,tslen,VR,VB,mu,Vini,fireini,sti,ptype,VT)
%NEURONS This program creates a time-series of a simple hypothesis for coupled neurons based on Galves-Locherbach stochastic model.
%
%input:
%- Synaptic weight matrix: w
%- Time-series length: tslen
%- Neuron definitions: VR (Reset voltage after firing),VB (baseline voltage),mu (current lakeage)
%- Initial neuron state: Vini,fireini
%- External stimulii: sti(1:nonodes,1:tslen)
%- Probability distribution type: ptype: 
%                                  'erf' for error function
%                                  'lif' for step function, which produces LIF neurons
%-Threshold voltage: VT (if ptype='erf' it's the zero of err function)
%output:
%- Neuron firing timeseries: neuronfire
%- Membrane potential result: V
%- Probability of a neuron firing: P
%
%example (two neurons, 1->2):
% Vini(1,1:2)=rand(1,2)-0.5;fireini(Vini>0)=1;fireini(Vini<=0)=0;sti(1:2,1:1e5)=0;
% [neuronfire,V,P]=neuron([0 0;1 0],1e5,-1,-0.2,0.1,Vini,fireini,sti,'erf',0);
%
%(C)Arthur Valencio 4 July 2018
%Update 11 November 2018

    nonodes=length(w(1,:));
    
    %set all neurons equal or distinct
    if length(VR)==1
        VR(1:nonodes)=VR;
    end
    if length(VB)==1
        VB(1:nonodes)=VB;
    end
    if length(mu)==1
        mu(1:nonodes)=mu;
    end
    
    %initial condition
    neuronfire(1:nonodes,1)=fireini(1:nonodes);
    V(1:nonodes,1)=Vini(1:nonodes);
    P(1:nonodes,1:tslen)=0;
    
    %calculate time-series
	for t=2:tslen	
		for node=1:nonodes
			if neuronfire(node,t-1)==1 %neuron-fire: reset to VR
				neuronfire(node,t)=0;
				V(node,t)=VR(node);
            else
					sumterm=0;
					for j=1:nonodes %sumterm gives the effect of the firing neighbour neurons
						sumterm=sumterm+w(node,j)*neuronfire(j,t-1);
                    end
                    %calculation of the potential:
					V(node,t)=mu(node)*(V(node,t-1)-VB(node))+VB(node)+sti(node,t)+sumterm;
                    %probability of a neuron firing:
                    if strcmp(ptype,'erf') %GL using error function
                        P(node,t)=calcprobperc(V(node,t),VT);
                    else %Leaky Integrate-and-Fire
                        if V(node,t)>=VT
                            P(node,t)=1;
                        else
                            P(node,t)=0;
                        end
                    end
                    %make neuron fire:
                    neuronfire(node,t)=coinflip(P(node,t));
            end
        end
    end
end

function [phi]= calcprobperc(V,VT) %calculates the firing probability function (erf)
	phi=0.5+0.5*erf(V-VT);
end

function x= coinflip(p) %defines if the neuron firing occurs given a probability p
	r=rand();
	if r<=p
		x=1;
    else
		x=0;
    end
end