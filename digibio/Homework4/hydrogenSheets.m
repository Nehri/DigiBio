 %protein = pdbread('3bvu.pdb', 'ModelNum', 1);
% 
% 
% 
% %% 
% clc;
% numHelices = length(protein.Helix);
% 
% counts = [0 0 0];
% 
% %for i = 1:numHelices
% for i = 1:2   
%     startLocation = protein.Helix(i).initSeqNum;
%     helixLength = protein.Helix(i).length;
%     endLocation = startLocation + helixLength;
%     
%     for j = startLocation:endLocation
%         
%         donorAcid = find([protein.Model.Atom.resSeq] == j);
%         donorArray = findstr('N', [protein.Model.Atom(donorAcid).AtomName]);
%         
%         donorIndex = donorArray(1);
%         donorAtomNumber = donorAcid(donorIndex);
%         donorAtom = protein.Model.Atom(donorAtomNumber);
%         vecD = [donorAtom.X donorAtom.Y donorAtom.Z];
%         
%         for k = [3:5]
%             if j+k > endLocation
%                 break
%             end
%            
%             acceptorAcid = find([protein.Model.Atom.resSeq] == j+k);
%             acceptorArray = findstr('O', [protein.Model.Atom(acceptorAcid).AtomName]);
%             
%             acceptorIndex = acceptorArray(1);
%             acceptorAtomNumber = acceptorAcid(acceptorIndex);
%             display(acceptorAtomNumber);
%             acceptorAtom = protein.Model.Atom(acceptorAtomNumber);
%             vecA = [acceptorAtom.X acceptorAtom.Y acceptorAtom.Z];
%             
%             
%             antecedentAtomNumber = acceptorAtomNumber - 1;
%             display(antecedentAtomNumber);
%             antecedentAtom = protein.Model.Atom(antecedentAtomNumber);
%             vecB = [antecedentAtom.X antecedentAtom.Y antecedentAtom.Z];
%             
%             hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k);
%             hydrogenArray = findstr('H', [protein.Model.Atom(hydrogenAcid).AtomName]);
%             
%             % Checks back one amino acid more in case no hydrogen was found
%             if isempty(hydrogenArray)
%                 hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k-1);
%                 hydrogenArray = findstr('H', [protein.Model.Atom(hydrogenAcid).AtomName]);
%             end
%             
%             hydrogenIndex = hydrogenArray(1);
%             hydrogenAtomNumber = hydrogenAcid(hydrogenIndex);
%             display(hydrogenAtomNumber);
%             hydrogenAtom = protein.Model.Atom(hydrogenAtomNumber);
%             vecH = [hydrogenAtom.X hydrogenAtom.Y hydrogenAtom.Z];
%             
%             result = hydrogen_analysis(vecD,vecA,vecH,vecB);
%             
%             if result
%                 counts(k-2) = counts(k-2)+1;
%             end
%         end
%     end 
% end
% 
% display(counts(1));
% display(counts(2));
% display(counts(3));


protein = pdbread('7AHL_B.pdb', 'ModelNum',1);


%%
clc;
numHelices = length(protein.Sheet);
 
counts = [0 0 0];
%for i = 1:numHelices
for i = 1:1
    startLocation = protein.Helix(i).initSeqNum;
    helixLength = protein.Helix(i).length;
    endLocation = startLocation + helixLength-1;
    
    for j = startLocation:endLocation
        donorAcid = find([protein.Model.Atom.resSeq] == j);
        donorArray = find([protein.Model.Atom(donorAcid).AtomName] == 'N');
        donorIndex = donorArray(1);
        donorAtom = donorAcid(donorIndex);
        donor = protein.Model.Atom(donorAtom);
        vecD = [donor.X donor.Y donor.Z];
        
        for k = 3:5
            if (j+k) > endLocation 
                break;
            end
            
            acceptorAcid = find([protein.Model.Atom.resSeq] == (j+k));
            acceptorArray = find([protein.Model.Atom(acceptorAcid).AtomName] == 'O');
            acceptorIndex = acceptorArray(1);
            acceptorAtom = acceptorAcid(acceptorIndex);
            acceptor = protein.Model.Atom(acceptorAtom);
            vecA = [acceptor.X acceptor.Y acceptor.Z];
            
            antecedentAtom = acceptorAtom - 1;
            antecedent = protein.Model.Atom(antecedentAtom);
            vecB = [antecedent.X antecedent.Y antecedent.Z];
            
            hydrogenAcid = find([protein.Model.Atom.resSeq] == (j+k));
            hydrogenArray = find([protein.Model.Atom(hydrogenAcid).AtomName] == 'H');
            
            % Checks back one amino acid more in case no hydrogen was found
            if isempty(hydrogenArray)
                hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k-1);
                hydrogenArray = find([protein.Model.Atom(hydrogenAcid).AtomName] == 'H');
            end
            
            hydrogenIndex = hydrogenArray(1);
            hydrogenAtom = hydrogenAcid(hydrogenIndex);
            hydro = protein.Model.Atom(hydrogenAtom);
            vecH = [hydro.X hydro.Y hydro.Z];
            
            result = hydrogen_analysis(vecD,vecA,vecH,vecB);
            
            if result
                display('hit');
                counts(k-2) = counts(k-2)+1;
            end

        end
    end
end