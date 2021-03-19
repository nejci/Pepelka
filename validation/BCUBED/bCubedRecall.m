function B3R = bCubedRecall(element, refrClust, respPart)
	respClust = (respPart == respPart(element)); % cluster of response partition that contains element
	nRefrClust = sum(refrClust);
	if  nRefrClust == 0
		B3R = 1;
	else
		B3R = sum(refrClust & respClust) / nRefrClust; 
	end
end