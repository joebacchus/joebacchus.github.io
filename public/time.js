const birth = new Date(2002, 3, 16, 0, 0, 0); 

function update() {
    const now = new Date();
    let diff = now - birth;

    const msperyear = 1000 * 60 * 60 * 24 * 365.25;
    const age = diff / msperyear;

    const years = Math.floor(age);
    const fractional = age - years;

    let remaining = fractional * msperyear;

    const msperday = 1000 * 60 * 60 * 24;
    const days = Math.floor(remaining / msperday);
    remaining -= days * msperday;

    const hours = Math.floor(remaining / (1000 * 60 * 60));
    remaining -= hours * 1000 * 60 * 60;

    const minutes = Math.floor(remaining / (1000 * 60));
    remaining -= minutes * 1000 * 60;

    const seconds = Math.floor(remaining / 1000);
    remaining -= seconds * 1000;

    const milliseconds = Math.floor(remaining); 

    const timeString = `${years}y ${String(days).padStart(3,'0')}d ${String(hours).padStart(2,'0')}h ${String(minutes).padStart(2,'0')}m ${String(seconds).padStart(2,'0')}s ${String(milliseconds).padStart(3,'0')}ms`;

    const clockElement = document.getElementById('clock');
    if (clockElement) {
        clockElement.textContent = timeString;
    }
}

update();
setInterval(update, 50);