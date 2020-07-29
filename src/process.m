clear all

function someval = read_plot(n,filename,yscale)
  dim=2+2*n;
  #filename='dfo_snm'
  j=0;
  for filno=1:30
    for d=1:2
      perturbstr=int2str(filno);
      if filno<10
        perturbstr=sprintf('%s%s',"0",perturbstr);
      end
      filenamelocal=sprintf('%s%s%s%s%d%s',filename, '_',perturbstr,'_',d);
      % Next section is a per-nucleus summary of the different Fayans computations
      % that were run to collect all the required theoretical results.
      fid = fopen([filenamelocal,'.txt']);
      if fid >= 0
        while 1
          line = fgetl(fid);
          if line==-1
             break
          end
          b=strsplit(line,",");
    
          if length(b)==2
            j=j+1;
            probstore(j,:)=cellfun(@str2num,b(1:end));
            k=0;
          else
            k=k+1;
            datastore(k,:,j)=cellfun(@str2num,b(1:dim));
          end
        end
      end
    end
  end

  for step=[0.1 0.2]
    prob=0;
    datastore(datastore==0)=NaN;
    figure(1); clf; hold on;
    #Loop over all available problems(up to 60)
    for j=1:size(probstore,1)
      if probstore(j,2)==step
        prob=prob+1;
        funval(:,prob)=datastore(:,2,j);
        best(1,prob) = datastore(1,2,j);
        grad(1,prob) = norm(datastore(1,n+3:dim,j));
        #Loop over all optimizer iterations in the problem
        for k=2:length(datastore(:,2,j))
            best(k,prob) = min(best(k-1,prob),datastore(k,2,j));
            grad(k,prob) = norm(datastore(k,n+3:dim,j));
        end
      end
    end
    figure(1)
    # Plot the function values from iteration 3 instead of 1
    # because the values of iteration 2 is very high in some
    # cases and leads to lower values not being plotted
    plot(funval(3:end,:))
    box on
    axis tight
    xlabel('iterations')
    ylabel('f value')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_1.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_1.png'])
    end

    figure(2)
    plot(best)
    box on
    axis tight
    set(gca,'yscale','linear')
    xlabel('iterations')
    ylabel('best f value')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_2.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_2.png'])
    end

    figure(3)
    plot(grad)
    box on
    axis tight
    set(gca,'yscale','log')
    xlabel('iterations')
    ylabel('gradient norm')
    title(sprintf('%s %f',filename, step))
    if step==0.1
      print('-dpng',['fig_',filename,'_','1','_3.png'])
    end
    if step==0.2
      print('-dpng',['fig_',filename,'_','2','_3.png'])
    end
  end
  some_val = 0; % could be scalar, vector, matrix etc.
end
read_plot(2,'bfg_snm','linear')
read_plot(2,'bfg_pnm','linear')
read_plot(4,'dfo_snm','linear')
read_plot(4,'dfo_pnm','linear')
