# for means
slide.mean = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# for median
slide.median = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = median(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# for standard deviation
slide.sd = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sd(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# for sum
slide.sum = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sum(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# for max
slide.max = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = max(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# for min
slide.min = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = min(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# fraction above
slide.count.gt = function(data, value, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sum(data[spots[i]:(spots[i]+window)]>value)
  }
  return(result)
}

# fraction below
slide.count.lt = function(data, value, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sum(data[spots[i]:(spots[i]+window)]<value)
  }
  return(result)
}

# fraction above or equal
slide.count.ge = function(data, value, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sum(data[spots[i]:(spots[i]+window)]>=value)
  }
  return(result)
}

# fraction below or equal
slide.count.le = function(data, value, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = sum(data[spots[i]:(spots[i]+window)]<=value)
  }
  return(result)
}

# concatenate values in a vector
slide.returnstring = function(data, sep, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] = paste(data[spots[i]:(spots[i]+window)], collapse=sep)
  }
  return(result)
}

# concatenate values in a vector
slide.returnvector = function(data, window, step){
  total = length(data)
  spots = seq(from=1, to=(total-window), by=step)
  result = list()
  for(i in 1:length(spots)){
    result[i] = c(data[data[spots[i]:(spots[i]+window)]])
  }
  return(result)
}
